#!/usr/bin/env python3
"""
pypomp_bootstrap_lrt.py
=======================
GPU / pypomp port of the Single-species **Dent SIRJPF** parametric-bootstrap-LRT
computation -- i.e. the GPU replacement for  fit_null.R  and  fit_alt.R  in this
directory.  One invocation fits ONE bootstrap dataset b under EITHER the null
(all-shared) OR one alternative (a single parameter block made unit-specific),
exactly as the two R scripts do, and writes a result consumed (after a one-time
R bridge -- see below) by the UNCHANGED  collect_lrt.R.

This is the bootstrap-LRT analogue of the validated SIRJPF2 ports
  ../../../../Mixed-species/SIRJPF2/model/pypomp_all_shared.py        (all-shared MPIF, GPU-confirmed: ll -880.44 vs R -880.56)
  ../../../../Mixed-species/SIRJPF2/coverage_study/pypomp_coverage_profile.py  (MPIF batching + R-free CSV I/O)
The JAX model below is the SIRJPF2 model REDUCED to single-species Dent (n-species
only); see the term-by-term mapping in the module docstring section "MODEL DERIVATION".

Usage
-----
    python pypomp_bootstrap_lrt.py <b 1..100> <target>

      target = "null"                       -> all-shared fit  (== fit_null.R)
      target in {xi, theta_Sn, theta_In,    -> make that one block unit-specific
                 theta_P, probn, rn, f_Sn}     and fit  (== fit_alt.R <b> <target>)

    # examples
    python pypomp_bootstrap_lrt.py 7 null
    python pypomp_bootstrap_lrt.py 7 theta_Sn

It reads simulated_data/sim_data_all.csv (the b-th panel) + true_params.csv with
pandas (NO R, NO rdata on the GPU node), fits via 2-round MPIF, pfilters, and
writes an .rds via pyreadr (matches the SIRJPF2 / Lum-SRJF reference ports):
  PYPOMP_OVERWRITE=0 (default): SIDECAR, does NOT touch the canonical R result
      null:  results_null/lrt_null_<b>_pypomp.rds
      alt:   results_alt/lrt_<target>_<b>_pypomp.rds
  PYPOMP_OVERWRITE=1:           canonical name collect_lrt.R reads directly
      null:  results_null/lrt_null_<b>.rds
      alt:   results_alt/lrt_<target>_<b>.rds
The .rds carries EXACTLY the fields collect_lrt.R reads (flattened to a one-row
data.frame): ll, se, coef.<name>, Np, Mp, Np_rep, block, Nmif.round1/round2.

-------------------------------  R-FREE I/O  ---------------------------------
* INPUT  (R-free): sim_data_all.csv + true_params.csv are produced ONCE from the
  committed .rds by  make_csv.R  (a tiny converter, embedded as a comment at the
  bottom of this file); the GPU node only reads CSVs.  Do NOT modify the .rds.
* OUTPUT: written as .rds via pyreadr.write_rds (one-row data.frame whose columns
  ARE collect_lrt.R's fields -- res$ll/res$Np/... index by name and return the
  length-1 value, so collect_lrt.R needs NO change).  PYPOMP_OVERWRITE gates the
  filename (sidecar vs canonical) so smoke runs never clobber the canonical R .rds.

================================  FIDELITY NOTE  ===============================
JAX's PRNG is NOT R's, so the particle draws / mif perturbations are not
bit-identical to fit_null.R / fit_alt.R; this reproduces the SAME estimator and
workflow (identical reduced-SIRJPF model, identical 2-round MPIF, Np=Mp=1500,
Np_rep=10, Nmif 150/250, cooling a=0.7, rw_sd, top-25%-replicate-x4 selection,
block=TRUE for alternatives).  Compare STATISTICALLY (the bootstrap Lambda
distribution and p_boot agree within Monte-Carlo error), not byte-for-byte.
R<->pypomp construct mapping annotated inline as  [R: fit_*.R <line/desc>].
===============================================================================
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time
from pathlib import Path

# ----------------------------------------------------------------------------
# 0. JAX device / precision configuration  -- MUST precede `import jax`
# ----------------------------------------------------------------------------
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")

_USE_GPU = os.environ.get("PYPOMP_USE_GPU", "1").lower() not in {"0", "false", "no", "cpu"}
if _USE_GPU:
    os.environ.pop("JAX_PLATFORMS", None)        # let JAX pick the GPU if present
else:
    os.environ["JAX_PLATFORMS"] = "cpu"

import jax  # noqa: E402

# float64: R/pomp likelihoods are double precision.  Keep x64 ON for fidelity.
_USE_X64 = os.environ.get("PYPOMP_X64", "1").lower() not in {"0", "false", "no"}
jax.config.update("jax_enable_x64", _USE_X64)

import jax.numpy as jnp  # noqa: E402
from jax.scipy.special import gammaln  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pypomp as pp  # noqa: E402


# ----------------------------------------------------------------------------
# 1. Constants  [R: fit_null.R / fit_alt.R]
# ----------------------------------------------------------------------------
N_UNITS = 8                                       # [R: for (i in 1:8); single-species Dent has 8 mesocosms]
UNIT_NAMES = [f"u{i}" for i in range(1, N_UNITS + 1)]

# Algorithmic params -- HARD-CODED to match fit_null.R / fit_alt.R exactly.
# [R: fit_null.R lines 124-126]   Np = 1500, Np_rep = 10, Mp = 1500
NP = 1500
NP_REP = 10
MP = 1500
NMIF_ROUND1 = 150                                 # [R: Nmif = 150]
NMIF_ROUND2 = 250                                 # [R: Nmif = 250]
RW_ROUND1 = 0.05                                  # [R: dent_rw.sd = 0.05]
RW_ROUND2 = 0.04                                  # [R: dent_rw.sd = 0.04]
COOLING_A = 0.7                                    # [R: cooling.fraction.50 = 0.7]
# Number of mif STARTS per round.  [R: 1:(3 * getDoParWorkers()) with 100 cores = 300].
# Each start is an independent fit (R's foreach), so this only sets multistart
# robustness; it can be tuned down for a smoke test via --n-starts.
N_STARTS_DEFAULT = 300

# A run-level smoke knob (does NOT affect the production run, which is run_level 3
# == the hard-coded NP/MP above).  run_level {1,2} shrink particles for CPU checks.
ALGORITHMIC_PARAMS = {
    "Np": [50, 500, NP],
    "Np_rep": [2, 10, NP_REP],
    "Mp": [50, 500, MP],
}

# Outer batch over mif STARTS (GPU-memory knob; batching does NOT change results --
# each start is independent, exactly like R's foreach).
BATCH_SIZE = int(os.environ.get("PYPOMP_BATCH", "150"))

# 15-parameter canonical order  [R: fit_null.R param_names line 90].
PARAM_NAMES = [
    "sigSn", "sigIn", "sigF", "sigP", "f_Sn", "rn", "k_In", "k_Sn",
    "sigJn", "probn", "theta_Sn", "theta_In", "theta_P", "theta_Jn", "xi",
]

# State vector  [R: fit_null.R state_names line 92].  Single-species Dent: n-species
# compartments only (no Si/Ii/Ji, no T_Si/T_Ii).
SIRJPF_STATENAMES = [
    "Sn", "In", "Jn", "F", "P", "T_Sn", "T_In", "error_count",
]

# Parameters carried on the LOG scale  [R: fit_null.R parameter_trans(log=...) lines 84-88].
# Every estimated parameter is log-transformed.  sigSn is fixed at 0 with rw.sd=0
# (it IS in the R log= list, but its value 0 is pinned by rw.sd=0; log(0)=-inf is
# avoided by keeping it out of the log transform here, exactly as the SIRJPF2 ports
# pin sigSn/sigSi at the boundary).
SIRJPF_LOG_PARAMS = (
    "sigIn", "sigF", "sigP", "f_Sn", "rn", "k_In", "k_Sn",
    "sigJn", "probn", "theta_Sn", "theta_In", "theta_P", "theta_Jn", "xi",
)
FIXED_ZERO_PARAMS = ("sigSn",)                    # [R: sigSn=0 with rw.sd sigSn=0]

# The seven alternatives this family tests  [R: fit_alt.R alt_params_map lines 18-26].
# Each maps to the single parameter block made unit-specific.
ALT_PARAMS_MAP = {
    "xi":       ["xi"],
    "theta_Sn": ["theta_Sn"],
    "theta_In": ["theta_In"],
    "theta_P":  ["theta_P"],
    "probn":    ["probn"],
    "rn":       ["rn"],
    "f_Sn":     ["f_Sn"],
}


# ============================================================================
#                            MODEL DERIVATION
# ============================================================================
# The JAX model below is the SIRJPF2 model (Mixed-species) REDUCED to the
# single-species Dent SIRJPF Csnippet in fit_null.R / fit_alt.R (lines 21-82),
# term-by-term:
#
# COMPARTMENTS DROPPED vs SIRJPF2 (which had Sn,In,Jn,Si,Ii,Ji,F,P):
#   * Si, Ii, Ji  (the invasive/i-species)  -> single-species Dent keeps only the
#     n-species Sn,In,Jn.  Accumvars T_Si, T_Ii dropped (no lum observables).
#   KEPT: Sn, In, Jn, F, P, T_Sn, T_In, error_count.
#
# PARAMETERS DROPPED vs SIRJPF2 (26 -> 15):
#   sigSi, sigIi, sigJi, theta_Si, theta_Ii, theta_Ji, f_Si, ri, probi, k_Si, k_Ii
#   (all the i-species params).  KEPT the 15 n-species params in PARAM_NAMES.
#
# RPROC, line-by-line  [R: fit_null.R lines 21-54 == fit_alt.R lines 55-88]:
#   noiSn = rnorm(0, sigSn*sqrt(dt))                      [R:26]  -> noiSn below
#   noiIn = rnorm(0, sigIn*sqrt(dt))                      [R:27]
#   noiJn = rnorm(0, sigJn*sqrt(dt))                      [R:28]
#   noiF  = rnorm(0, sigF *sqrt(dt))                      [R:29]
#   noiP  = rnorm(0, sigP *sqrt(dt))                      [R:30]
#   Sn_term = 0.1*Jn*dt - theta_Sn*Sn*dt - probn*f_Sn*Sn*P*dt - delta*Sn*dt + Sn*noiSn   [R:32]
#   Jn_term = rn*f_Sn*F*Sn*dt - 0.1*Jn*dt - theta_Jn*Jn*dt - delta*Jn*dt + Jn*noiJn       [R:33]
#   In_term = probn*f_Sn*Sn*P*dt - theta_In*In*dt - delta*In*dt + In*noiIn                [R:34]
#   F_term  = F*noiF - f_Sn*F*(Sn + xi*In + 1*Jn)*dt - delta*F*dt + 0.37*dt               [R:35]
#             (NOTE: NO i-species term in F_term -- the SIRJPF2 had a second -f_Si*F*(...) )
#   P_term  = 30*theta_In*In*dt - f_Sn*(Sn + xi*In)*P*dt - theta_P*P*dt - delta*P*dt + P*noiP  [R:36]
#             (NOTE: SIRJPF2 had +30*theta_Ii*Ii and -f_Si*(Si+xi*Ii)*P; both dropped.)
#   Euler update  F+=,Sn+=,In+=,Jn+=,P+=                  [R:38-42]
#   parasite pulse:  if |t-4|<0.001  P += 25              [R:44]
#   soft constraints (delta=0.013, lambda_J=0.1, mu_food=0.37 are HARD-CODED):
#     Sn<0|>1e5 -> Sn=0, error_count += 1                 [R:46]
#     F <0|>1e20-> F =0, error_count += 1000               [R:47]   (SIRJPF2 used 1e3)
#     In<0|>1e5 -> In=0, error_count += 0.001              [R:48]
#     Jn<0|>1e5 -> Jn=0, error_count += 0.001              [R:49]
#     (P<0|>1e20)&(t>3.9) -> P=0, error_count += 1e-6      [R:50]
#   T_Sn = fabs(Sn);  T_In = fabs(In)                      [R:52-53]
#   NOTE the R clamps the violating state to 0 (not clip-to-bounds); we replicate
#   "set to 0 on violation" via jnp.where, matching the C exactly.
#
# RINIT  [R: fit_null.R lines 56-65]:
#   Sn = 3   (NB: SIRJPF2 used Sn=2.333; Dent uses 3),  F = 16.667,  rest 0.
#
# DMEAS (negative binomial, 2 observables)  [R: fit_null.R lines 67-77]:
#   if error_count>0 -> lik = -150
#   else lik = dnbinom_mu(dentadult, k_Sn, T_Sn) + dnbinom_mu(dentinf, k_In, T_In)   (log scale)
#   (SIRJPF2 had 4 terms incl. lumadult/luminf; here only the 2 dent observables.)
#
# dt = 0.25  [R: euler(delta.t = 1/4)],  accumvars = ("error_count",)  [R: accumvars].
# ============================================================================


def sirjpf_rproc(X_, theta_, key, covars, t, dt):
    """One Euler-Maruyama step  [R: fit_null.R dyn_rpro lines 21-54]."""
    Sn, Jn, In = X_["Sn"], X_["Jn"], X_["In"]
    F, P = X_["F"], X_["P"]
    error_count = X_["error_count"]

    sigSn, sigIn = theta_["sigSn"], theta_["sigIn"]
    sigJn = theta_["sigJn"]
    sigF, sigP = theta_["sigF"], theta_["sigP"]
    theta_Sn, theta_In = theta_["theta_Sn"], theta_["theta_In"]
    theta_Jn = theta_["theta_Jn"]
    theta_P = theta_["theta_P"]
    f_Sn = theta_["f_Sn"]
    rn = theta_["rn"]
    probn = theta_["probn"]
    xi = theta_["xi"]

    delta = 0.013          # [R: double delta = 0.013]
    mu_food = 0.37         # [R: + 0.37 * dt]
    lambda_J = 0.1         # [R: 0.1 * Jn ... maturation rate]

    keys = jax.random.split(key, 5)
    sqdt = jnp.sqrt(dt)
    noiSn = sigSn * sqdt * jax.random.normal(keys[0])     # [R:26]
    noiIn = sigIn * sqdt * jax.random.normal(keys[1])     # [R:27]
    noiJn = sigJn * sqdt * jax.random.normal(keys[2])     # [R:28]
    noiF = sigF * sqdt * jax.random.normal(keys[3])       # [R:29]
    noiP = sigP * sqdt * jax.random.normal(keys[4])       # [R:30]

    # [R:32-36]  -- single-species (n) flows; F_term / P_term carry NO i-species term.
    Sn_term = (lambda_J * Jn * dt - theta_Sn * Sn * dt
               - probn * f_Sn * Sn * P * dt - delta * Sn * dt + Sn * noiSn)
    Jn_term = (rn * f_Sn * F * Sn * dt - lambda_J * Jn * dt
               - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn)
    In_term = (probn * f_Sn * Sn * P * dt - theta_In * In * dt
               - delta * In * dt + In * noiIn)
    F_term = (F * noiF - f_Sn * F * (Sn + xi * In + 1.0 * Jn) * dt
              - delta * F * dt + mu_food * dt)
    P_term = (30.0 * theta_In * In * dt
              - f_Sn * (Sn + xi * In) * P * dt
              - theta_P * P * dt - delta * P * dt + P * noiP)

    Sn_new = Sn + Sn_term     # [R:39]
    In_new = In + In_term     # [R:40]
    Jn_new = Jn + Jn_term     # [R:41]
    F_new = F + F_term        # [R:38]
    P_new = P + P_term        # [R:42]
    # [R:44]  parasite inoculation pulse at day 4
    P_new = P_new + jnp.where(jnp.abs(t - 4.0) < 0.001, 25.0, 0.0)

    # [R:46-50]  soft constraints -> error_count, then the VIOLATING state set to 0
    # (the C sets the state to 0, NOT clip-to-bound -- replicate exactly).
    def viol(x, hi):
        return (x < 0.0) | (x > hi)

    v_Sn = viol(Sn_new, 1e5)
    v_F = viol(F_new, 1e20)
    v_In = viol(In_new, 1e5)
    v_Jn = viol(Jn_new, 1e5)
    v_P = viol(P_new, 1e20) & (t > 3.9)

    eps = 0.0
    eps += jnp.where(v_Sn, 1.0, 0.0)          # [R:46]  error_count += 1
    eps += jnp.where(v_F, 1000.0, 0.0)        # [R:47]  error_count += 1000
    eps += jnp.where(v_In, 0.001, 0.0)        # [R:48]  error_count += 0.001
    eps += jnp.where(v_Jn, 0.001, 0.0)        # [R:49]  error_count += 0.001
    eps += jnp.where(v_P, 1.0e-6, 0.0)        # [R:50]  error_count += 0.000001

    Sn_new = jnp.where(v_Sn, 0.0, Sn_new)
    F_new = jnp.where(v_F, 0.0, F_new)
    In_new = jnp.where(v_In, 0.0, In_new)
    Jn_new = jnp.where(v_Jn, 0.0, Jn_new)
    P_new = jnp.where(v_P, 0.0, P_new)

    return {
        "Sn": Sn_new, "In": In_new, "Jn": Jn_new, "F": F_new, "P": P_new,
        "T_Sn": jnp.abs(Sn_new), "T_In": jnp.abs(In_new),   # [R:52-53]
        "error_count": error_count + eps,                   # accumvar -> reset each obs
    }


def sirjpf_rinit(theta_, key, covars, t0):
    """Deterministic initial state  [R: fit_null.R dyn_init lines 56-65]."""
    return {
        "Sn": jnp.array(3.0), "In": jnp.array(0.0), "Jn": jnp.array(0.0),
        "F": jnp.array(16.667), "P": jnp.array(0.0),
        "T_Sn": jnp.array(0.0), "T_In": jnp.array(0.0),
        "error_count": jnp.array(0.0),
    }


def _nb_logpmf(y, mu, size):
    """Negative-binomial log-pmf in (mu, size) form  [R: dnbinom_mu(y, size, mu)]."""
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    return (gammaln(y + size) - gammaln(size) - gammaln(y + 1.0)
            + size * jnp.log(size / (size + mu))
            + y * jnp.log(mu / (size + mu)))


def sirjpf_dmeas(Y_, X_, theta_, covars, t):
    """2-D NB measurement log-likelihood  [R: fit_null.R dmeas lines 67-77]."""
    ll = (_nb_logpmf(Y_["dentadult"], X_["T_Sn"], theta_["k_Sn"])
          + _nb_logpmf(Y_["dentinf"], X_["T_In"], theta_["k_In"]))
    # [R: if (error_count > 0.0) lik = -150]
    return jnp.where(X_["error_count"] > 0.0, -150.0, ll)


def _nb_sample(key, mu, size):
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    k1, k2 = jax.random.split(key)
    g = jax.random.gamma(k1, size) * (mu / size)   # Gamma-Poisson = NB
    return jax.random.poisson(k2, g)


def sirjpf_rmeas(X_, theta_, key, covars, t):
    """Measurement simulator  [R: fit_null.R rmeas lines 79-82]  (unused by LRT fit)."""
    keys = jax.random.split(key, 2)
    return jnp.array([
        _nb_sample(keys[0], X_["T_Sn"], theta_["k_Sn"]),
        _nb_sample(keys[1], X_["T_In"], theta_["k_In"]),
    ], dtype=float)


def sirjpf_to_est(theta):
    """natural -> estimation scale  [R: parameter_trans(log=...) toEst]."""
    out = {**theta}
    for name in SIRJPF_LOG_PARAMS:
        out[name] = jnp.log(jnp.maximum(theta[name], 1e-30))
    for name in FIXED_ZERO_PARAMS:
        out[name] = theta[name]        # pinned at 0, not log-transformed
    return out


def sirjpf_from_est(theta):
    """estimation -> natural scale  [R: parameter_trans fromEst]."""
    out = {**theta}
    for name in SIRJPF_LOG_PARAMS:
        out[name] = jnp.exp(theta[name])
    for name in FIXED_ZERO_PARAMS:
        out[name] = theta[name]
    return out


sirjpf_par_trans = pp.ParTrans(to_est=sirjpf_to_est, from_est=sirjpf_from_est)


# ----------------------------------------------------------------------------
# 2. R-FREE I/O: read pre-converted CSVs with pandas; write result .rds (pyreadr).
#    sim_data_all.csv (b,unit,day,dentadult,dentinf) + true_params.csv (name,value)
#    are produced ONCE from the .rds by make_csv.R (embedded at the bottom).
# ----------------------------------------------------------------------------
def read_sim_data(sim_dir: Path, b: int) -> dict[str, pd.DataFrame]:
    """Read the b-th simulated panel from sim_data_all.csv (cols: b,unit,day,2 obs)."""
    df = pd.read_csv(sim_dir / "sim_data_all.csv")
    df = df[df["b"] == b]
    if df.empty:
        raise ValueError(f"no rows for dataset b={b} in {sim_dir/'sim_data_all.csv'}")
    out = {}
    for unit_name in UNIT_NAMES:
        out[unit_name] = (df[df["unit"] == unit_name]
                          [["day", "dentadult", "dentinf"]]
                          .astype(float).sort_values("day").reset_index(drop=True))
    return out


def read_params(sim_dir: Path) -> dict[str, float]:
    """Read true_params.csv (name,value)."""
    tp = pd.read_csv(sim_dir / "true_params.csv")
    params = dict(zip(tp["name"].astype(str), tp["value"].astype(float)))
    missing = [n for n in PARAM_NAMES if n not in params]
    if missing:
        raise ValueError(f"true_params.csv missing parameters: {missing}")
    return params


def write_result_rds(path: Path, ll: float, se: float, coef: dict[str, float],
                     Np: int, Mp: int, Np_rep: int) -> None:
    """Write the result as an .rds in collect_lrt.R's schema (matches SIRJPF2/Lum ports).

    collect_lrt.R reads res$ll, res$se, and res$Np/Mp/Np_rep/block (metadata guard).
    fit_null.R/fit_alt.R save list(ll, se, coef, Np, Mp, Np_rep, block,
    Nmif=c(round1=150, round2=250)).  pyreadr.write_rds cannot write an arbitrary R
    named-list, so we write a ONE-ROW data.frame whose columns ARE those names:
    collect_lrt.R's `res$ll`, `res$Np`, ... index by name, which on a 1-row
    data.frame returns the length-1 value -> identical numerics, NO collect_lrt.R
    change needed.  coef is flattened into coef.<name> columns; Nmif into
    Nmif.round1 / Nmif.round2."""
    import pyreadr
    path.parent.mkdir(parents=True, exist_ok=True)
    rec = {
        "ll": [float(ll)],
        "se": [float(se)],
        "Np": [int(Np)],
        "Mp": [int(Mp)],
        "Np_rep": [int(Np_rep)],
        "block": [True],
        "Nmif.round1": [int(NMIF_ROUND1)],
        "Nmif.round2": [int(NMIF_ROUND2)],
    }
    for k, v in coef.items():
        rec[f"coef.{k}"] = [float(v)]
    pyreadr.write_rds(str(path), pd.DataFrame(rec))


# ----------------------------------------------------------------------------
# 3. rw_sd  [R: fit_null.R / fit_alt.R rw_sd(...) lines 142-146 / 182-186]
#    ALL params perturbed at `value` EXCEPT sigSn = 0.  (The alt block is also
#    perturbed at `value`; block=TRUE makes those per-unit moves.)
# ----------------------------------------------------------------------------
def make_rw_sd(value: float) -> dict[str, float]:
    rw = {name: value for name in PARAM_NAMES}
    rw["sigSn"] = 0.0          # [R: sigSn = 0]
    return rw


# ----------------------------------------------------------------------------
# 4. Panel + one MPIF round (batched over starts for GPU memory)
# ----------------------------------------------------------------------------
def logmeanexp_se(values: np.ndarray, axis: int):
    """logmeanexp and its Monte-Carlo SE along `axis`  [matches pomp::logmeanexp se=TRUE].
    [R: panel_logmeanexp(ll, MARGIN=1, se=TRUE)]."""
    m = np.nanmax(values, axis=axis, keepdims=True)
    w = np.exp(values - m)
    n = np.sum(~np.isnan(values), axis=axis)
    mean_w = np.nanmean(w, axis=axis)
    est = np.squeeze(m, axis=axis) + np.log(mean_w)
    se = np.nanstd(w, axis=axis, ddof=1) / mean_w / np.sqrt(n)
    return est, se


def build_pomp_dict(sim_data: dict[str, pd.DataFrame], theta: dict[str, float]):
    """8 per-unit pp.Pomp objects  [R: fit_null.R pomplist loop lines 97-119]."""
    pomp_dict = {}
    for unit_name in UNIT_NAMES:
        ys_u = (sim_data[unit_name]
                .set_index("day")[["dentadult", "dentinf"]]
                .astype(float))
        pomp_dict[unit_name] = pp.Pomp(
            ys=ys_u, theta=pp.PompParameters(theta), statenames=SIRJPF_STATENAMES, t0=1.0,
            rinit=sirjpf_rinit, rproc=sirjpf_rproc,
            dmeas=sirjpf_dmeas, rmeas=sirjpf_rmeas,
            par_trans=sirjpf_par_trans, dt=0.25,
            accumvars=("error_count",),           # [R: accumvars = c("error_count")]
        )
    return pomp_dict


def panel_theta_to_row(theta_dict: dict, specific_names: list[str]) -> dict[str, float]:
    """
    Extract a flat param dict from one panel-theta spec, for the NEXT round's start
    table.  Shared params -> bare name.  Unit-specific params -> per-unit columns
    "<p>[u1]".."<p>[u8]".  (For null, specific_names is empty -> 15 bare columns.)
    """
    shared_df = theta_dict["shared"]
    out = {n: float(shared_df.loc[n, "shared"]) for n in shared_df.index}
    if specific_names:
        unit_df = theta_dict["unit_specific"]
        for p in specific_names:
            for u in UNIT_NAMES:
                out[f"{p}[{u}]"] = float(unit_df.loc[p, u])
    return out


def panel_theta_to_coef(theta_dict: dict, specific_names: list[str]) -> dict[str, float]:
    """
    Build the panelPomp-style coef() naming for the OUTPUT JSON:
      shared params -> bare name;  unit-specific -> "<p>[u1]".."<p>[u8]".
    Matches the R results' coef names exactly (verified against an existing
    results_alt/lrt_theta_Sn_*.rds: bare shared names + "theta_Sn[u1]"..[u8]).
    """
    return panel_theta_to_row(theta_dict, specific_names)


def _run_mif_batch(pomp_dict, panel_params, specific_names, rw_sd_dict, nmif,
                   mif_particles, pfilter_particles, pfilter_reps,
                   mif_key, pf_key, vmap_chunk_size):
    """MPIF + replicated pfilter on ONE batch of starts -> list of (row, ll, se)."""
    panel = pp.PanelPomp(Pomp_dict=pomp_dict, theta=panel_params)

    # [R: mif2(Nmif, rw.sd, cooling.fraction.50=0.7, Np=Mp, block=TRUE)]
    # block=True is pypomp's MPIF: with a non-empty unit_specific (the ALT) it
    # block-resamples the per-unit columns; with an empty unit_specific (the NULL)
    # it is a no-op == PIF, matching the R "block=TRUE (no-op for all-shared null)".
    # pypomp >=0.4.6: cooling via RWSigma.geometric_cooling(a=...); a=0.7 == R 0.7.
    panel.mif(
        J=mif_particles, M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
        key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size,
    )

    # [R: ll <- replicate(Np_rep, unitLogLik(pfilter(m1, Np=Np)));
    #     panel_logmeanexp(ll, MARGIN=1, se=TRUE)]
    panel.pfilter(
        J=pfilter_particles, reps=pfilter_reps, key=pf_key,
        chunk_size=vmap_chunk_size if vmap_chunk_size else 1,
    )

    loglik_array = np.asarray(panel.results_history[-1].logLiks.values)  # (Nstarts,U,reps)
    assert loglik_array.ndim == 3, f"expected (Nstarts,U,reps), got {loglik_array.shape}"
    # logmeanexp over reps PER UNIT, THEN sum over units (these do not commute).
    unit_est, unit_se = logmeanexp_se(loglik_array, axis=2)   # (Nstarts, U) each
    panel_ll = np.nansum(unit_est, axis=1)                    # (Nstarts,)
    panel_se = np.sqrt(np.nansum(unit_se ** 2, axis=1))       # combine unit SEs

    results = []
    for idx, theta_dict in enumerate(panel.theta._to_list()):
        row = panel_theta_to_row(theta_dict, specific_names)
        results.append((row, float(panel_ll[idx]), float(panel_se[idx])))
    return results


def run_mif_round(pomp_dict, panel_params, specific_names, rw_sd_dict, nmif,
                  mif_particles, pfilter_particles, pfilter_reps, key_seed,
                  vmap_chunk_size):
    """
    One round over all starts, fed to pypomp in batches of BATCH_SIZE.
    Each start is an independent fit (R's foreach), so batching only caps peak GPU
    memory; results match an all-at-once run up to Monte-Carlo noise.  `panel_params`
    is a list of per-start theta specs (one PanelParameters built per batch).
    Returns a list of (row_dict, ll, se).
    """
    spec_list = panel_params  # list of {shared, unit_specific} dicts
    n = len(spec_list)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    results = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch_specs = spec_list[lo:lo + BATCH_SIZE]
        print(f"    batch {bi + 1}/{n_batches}  (starts {lo}..{lo + len(batch_specs) - 1})",
              flush=True)
        results.extend(_run_mif_batch(
            pomp_dict, pp.PanelParameters(theta=batch_specs), specific_names,
            rw_sd_dict, nmif, mif_particles, pfilter_particles, pfilter_reps,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size,
        ))
    return results


def select_and_replicate(results: list, specific_names: list[str]) -> list[dict]:
    """Top ceil(N/4) by ll, each replicated 4x  [R: fit_*.R select + rep(...,each=4)]."""
    ordered = sorted(results, key=lambda r: r[1], reverse=True)
    n_select = math.ceil(len(results) / 4.0)
    selected = ordered[:n_select]
    specs = []
    for row, _ll, _se in selected:
        for _ in range(4):
            specs.append(_row_to_spec(row, specific_names))
    return specs


def _row_to_spec(row: dict, specific_names: list[str]) -> dict:
    """Flat (possibly per-unit) row dict -> one {shared, unit_specific} theta spec."""
    shared_param_names = [n for n in PARAM_NAMES if n not in specific_names]
    shared_df = pd.DataFrame(
        {"shared": [float(row[name]) for name in shared_param_names]},
        index=shared_param_names,
    )
    if specific_names:
        unit_df = pd.DataFrame(
            {u: [float(row[f"{p}[{u}]"]) for p in specific_names] for u in UNIT_NAMES},
            index=specific_names,
        )
    else:
        unit_df = pd.DataFrame(index=[], columns=UNIT_NAMES)
    return {"shared": shared_df, "unit_specific": unit_df}


def initial_specs(theta0: dict[str, float], specific_names: list[str],
                  n_starts: int) -> list[dict]:
    """Round-1 starts: all = the warm start theta0  [R: shared.start=true_params,
    specific.start=rep(true_params[p],8)].  n_starts identical specs."""
    shared_param_names = [n for n in PARAM_NAMES if n not in specific_names]
    shared_df = pd.DataFrame(
        {"shared": [float(theta0[name]) for name in shared_param_names]},
        index=shared_param_names,
    )
    if specific_names:
        unit_df = pd.DataFrame(
            {u: [float(theta0[p]) for p in specific_names] for u in UNIT_NAMES},
            index=specific_names,
        )
    else:
        unit_df = pd.DataFrame(index=[], columns=UNIT_NAMES)
    spec = {"shared": shared_df, "unit_specific": unit_df}
    return [dict(spec) for _ in range(n_starts)]


# ----------------------------------------------------------------------------
# 5. Main: 2-round MPIF (== fit_null.R / fit_alt.R)
# ----------------------------------------------------------------------------
def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="GPU/pypomp bootstrap-LRT fit for ONE (dataset b, target) of "
                    "the single-species Dent SIRJPF family.")
    p.add_argument("b", type=int, help="bootstrap dataset index 1..100")
    p.add_argument("target", type=str,
                   help='"null" (all-shared) or one alternative: '
                        + ", ".join(ALT_PARAMS_MAP))
    p.add_argument("--run-level", type=int,
                   default=int(os.environ.get("PYPOMP_RUN_LEVEL", "3")),
                   choices=[1, 2, 3],
                   help="particle/iteration budget; 3 == production (Np=Mp=1500)")
    p.add_argument("--n-starts", type=int,
                   default=int(os.environ["PYPOMP_NSTARTS"]) if os.environ.get("PYPOMP_NSTARTS")
                   else N_STARTS_DEFAULT,
                   help="mif starts per round (R uses 3*100=300)")
    return p.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    b = args.b
    target = args.target
    rl = args.run_level
    if not (1 <= b <= 100):
        raise ValueError("dataset index must be 1..100")
    is_null = (target == "null")
    if not is_null and target not in ALT_PARAMS_MAP:
        raise ValueError(f"unknown target {target!r}; valid: null, "
                         + ", ".join(ALT_PARAMS_MAP))
    specific_names = [] if is_null else ALT_PARAMS_MAP[target]

    rl_i = rl - 1
    Mp = ALGORITHMIC_PARAMS["Mp"][rl_i]
    Np = ALGORITHMIC_PARAMS["Np"][rl_i]
    Np_rep = ALGORITHMIC_PARAMS["Np_rep"][rl_i]
    vmap_chunk = os.environ.get("PYPOMP_VMAP_CHUNK")
    vmap_chunk = int(vmap_chunk) if vmap_chunk else None

    print(f"=== pypomp Dent-SIRJPF bootstrap-LRT: b={b} target={target} "
          f"(run_level {rl}) ===")
    print(f"JAX backend = {jax.default_backend()}  x64 = {_USE_X64}  "
          f"devices = {jax.devices()}")
    print(f"Mp={Mp} Np={Np} Np_rep={Np_rep} n_starts={args.n_starts} "
          f"batch={BATCH_SIZE}  unit_specific={specific_names or '[] (all-shared null)'}")

    base = Path(__file__).resolve().parent
    sim_dir = base / "simulated_data"
    # PYPOMP_OVERWRITE=0 (default) -> SIDECAR "*_pypomp.rds" (does NOT touch the
    # canonical R lrt_*.rds); =1 -> write the canonical lrt_*.rds feeding collect_lrt.R.
    overwrite = os.environ.get("PYPOMP_OVERWRITE", "0").lower() in {"1", "true", "yes"}
    if is_null:
        name = f"lrt_null_{b}" if overwrite else f"lrt_null_{b}_pypomp"
        out_path = base / "results_null" / f"{name}.rds"
    else:
        name = f"lrt_{target}_{b}" if overwrite else f"lrt_{target}_{b}_pypomp"
        out_path = base / "results_alt" / f"{name}.rds"

    sim_data = read_sim_data(sim_dir, b)
    true_params = read_params(sim_dir)
    theta0 = {n: true_params[n] for n in PARAM_NAMES}
    pomp_dict = build_pomp_dict(sim_data, theta0)

    # ---- Round 1 ----  [R: Nmif=150, dent_rw.sd=0.05, warm start]
    print("Round 1 starting...", flush=True)
    _t0 = time.time()
    round1_specs = initial_specs(theta0, specific_names, args.n_starts)
    round1 = run_mif_round(
        pomp_dict, round1_specs, specific_names, make_rw_sd(RW_ROUND1),
        nmif=NMIF_ROUND1, mif_particles=Mp, pfilter_particles=Np,
        pfilter_reps=Np_rep, key_seed=801_000 + b, vmap_chunk_size=vmap_chunk,
    )
    print(f"Round 1 done. ({(time.time() - _t0) / 60:.1f} min)  "
          f"best ll = {max(r[1] for r in round1):.3f}", flush=True)

    # ---- Select top 25%, replicate x4 ----  [R: order/select + rep(...,each=4)]
    round2_specs = select_and_replicate(round1, specific_names)

    # ---- Round 2 ----  [R: Nmif=250, dent_rw.sd=0.04]
    print("Round 2 starting...", flush=True)
    _t1 = time.time()
    round2 = run_mif_round(
        pomp_dict, round2_specs, specific_names, make_rw_sd(RW_ROUND2),
        nmif=NMIF_ROUND2, mif_particles=Mp, pfilter_particles=Np,
        pfilter_reps=Np_rep, key_seed=801_000 + b + 100_000, vmap_chunk_size=vmap_chunk,
    )
    print(f"Round 2 done. ({(time.time() - _t1) / 60:.1f} min)  "
          f"[job total {(time.time() - _t0) / 60:.1f} min]", flush=True)

    # ---- Best of round 2  [R: which.max(lls[1,]); ll/se/coef] ----
    best_row, best_ll, best_se = max(round2, key=lambda r: r[1])
    coef = panel_theta_to_coef(_row_to_spec(best_row, specific_names), specific_names)

    write_result_rds(out_path, best_ll, best_se, coef, Np, Mp, Np_rep)
    print(f"\nSaved {out_path}")
    print(f"  ll = {best_ll:.4f}   se = {best_se:.4f}   (coef: {len(coef)} entries; "
          f"shared {len([k for k in coef if '[' not in k])} + "
          f"unit-specific {len([k for k in coef if '[' in k])})")
    if not overwrite:
        print("  NOTE: SIDECAR write (PYPOMP_OVERWRITE=0); canonical lrt_*.rds untouched. "
              "Set PYPOMP_OVERWRITE=1 to write the canonical name collect_lrt.R reads.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))


# ============================================================================
#  EMBEDDED ONE-TIME R BRIDGE (run LOCALLY where R is available; NOT on the
#  R-free GPU node).  Copy to its own .R file, or run via Rscript -e.
# ============================================================================
#
# ---- make_csv.R : ONE-TIME input converter (.rds -> .csv) ------------------
#   Produces simulated_data/sim_data_all.csv and simulated_data/true_params.csv,
#   which this script reads R-free on the GPU node.  (Already run for this dir.)
#
#   sim_dir <- "simulated_data"
#   rds_files <- sort(Sys.glob(file.path(sim_dir, "sim_data_*.rds")))
#   all_rows <- list()
#   for (f in rds_files) {
#     b <- as.integer(sub(".*sim_data_0*([0-9]+)\\.rds$", "\\1", f))
#     sd <- readRDS(f)
#     for (u in names(sd)) {
#       d <- sd[[u]]; colnames(d) <- c("day","dentadult","dentinf")
#       d$b <- b; d$unit <- u
#       all_rows[[length(all_rows)+1]] <- d[, c("b","unit","day","dentadult","dentinf")]
#     }
#   }
#   write.csv(do.call(rbind, all_rows), file.path(sim_dir,"sim_data_all.csv"), row.names=FALSE)
#   tp <- readRDS(file.path(sim_dir,"true_params.rds"))
#   write.csv(data.frame(name=names(tp), value=as.numeric(tp)),
#             file.path(sim_dir,"true_params.csv"), row.names=FALSE)
#
# Output is written DIRECTLY as .rds via pyreadr (see write_result_rds); no
# json->rds bridge is needed.  collect_lrt.R reads res$ll/res$se/res$Np/... by
# name off the one-row data.frame UNCHANGED.
# ============================================================================
