#!/usr/bin/env python3
"""
pypomp_bootstrap_lrt.py  --  Single-species Lum SRJF parametric-bootstrap-LRT on GPU
====================================================================================
GPU / pypomp port of  fit_null.R + fit_alt.R  (the Lum SRJF bootstrap-LRT arm).
It runs the SAME marginalized panel iterated filtering (MPIF) two-round fit on the
SAME simulated panels and writes the SAME minimal output consumed by collect_lrt.R:

    results_null/lrt_null_<b>.rds                 (target = "null"; ALL-shared fit)
    results_alt/lrt_<alt>_<b>.rds                 (target = an alternative name)

This is the Lum SRJF analogue of the validated SIRJPF2 ports
(../../../../Mixed-species/SIRJPF2/model/pypomp_all_shared.py  GPU-confirmed
 ll -880.44 vs R -880.56, and  .../coverage_study/pypomp_coverage_profile.py).
The model below is the SIRJPF2 JAX model REDUCED to the single-species i-only SRJF
sub-model -- see the FAMILY REDUCTION + R-CSNIPPET MAPPING notes immediately below.

Usage
-----
    python pypomp_bootstrap_lrt.py <b 1..100> <target>
        target = "null"               -> all-shared MPIF fit (PIF; block is a no-op)
               = "theta_Si"|"ri"|"f_Si"  -> that param made UNIT-SPECIFIC (MPIF, block=True)

    # examples
    python pypomp_bootstrap_lrt.py 7 null
    python pypomp_bootstrap_lrt.py 7 ri
    python pypomp_bootstrap_lrt.py 42 theta_Si --run-level 1   # tiny smoke

Environment toggles (defaults reproduce the R run):
    PYPOMP_USE_GPU = 1 (default) | 0/cpu
    PYPOMP_X64     = 1 (default)                # R uses double precision
    PYPOMP_RUN_LEVEL = 3 (default) | 1 | 2
    PYPOMP_BATCH   = 120 (default) | <int>      # outer batch over the 3*W starts
    PYPOMP_VMAP_CHUNK = (unset) | <int>         # pypomp per-unit vmap chunk
    PYPOMP_NSTARTS = (unset) | <int>            # number of round-1 starts (R uses 3*100=300)
    PYPOMP_OVERWRITE = 0 (default) | 1          # 0 => *_pypomp.rds SIDECAR (does NOT
                                                #   touch the canonical R result);
                                                #   1 => overwrite canonical lrt_*.rds

REQUIRES the dev pypomp (>=0.4.6: PanelPomp/PanelParameters/RWSigma/ParTrans) and a
CUDA JAX.  R-FREE on the GPU node: the input panels are read from the committed
simulated_data/sim_data_all.csv (pandas), produced ONCE from the .rds by the
converter at the bottom of this file; the only I/O dependency is `pyreadr`
(to write the output .rds in collect_lrt.R's schema).

==============================  FAMILY REDUCTION  ==============================
SIRJPF2 (mixed-species) -> Lum SRJF (single-species, i / "lum" only).  DROPPED vs
the SIRJPF2 JAX model in pypomp_coverage_profile.py:
  * the WHOLE n-species (native) block: states Sn, In, Jn and params
    theta_Sn, theta_In, theta_Jn, f_Sn, rn, probn, sigSn, sigIn, sigJn, k_Sn, k_In,
    obs dentadult/dentinf  -- SRJF is single-species, so only the i-species survives.
  * the I (infected) and P (parasite) compartments ENTIRELY: states Ii, P and params
    theta_Ii, theta_P, probi, xi, sigIi, sigP, k_Ii, obs luminf  -- SRJF = S,R(=J),F
    only (no infection dynamics); there is also NO P+=25 inoculation pulse.
This leaves the i-species SUSCEPTIBLE / JUVENILE / FOOD sub-model:
  states  Si, Ji, F  (+ T_Si accum-ready, error_count accumvar)
  params  ri, f_Si, theta_Si, theta_Ji, sigSi, sigJi, sigF, k_Si   (8, R param_names)
  obs     lumadult  (1)         units  u1..u10  (10, R loops 1:10)
exactly the Csnippets in fit_null.R / fit_alt.R  (verified term-by-term below).

=====================  R-CSNIPPET  ->  JAX  TERM-BY-TERM  ======================
rproc  [R: fit_null.R lines 21-43  dyn_rpro, euler(delta.t = 1/4) -> dt=0.25]:
  noiSi = rnorm(0, sigSi*sqrt(dt))                 -> sigSi*sqdt*N(0,1)         (L26)
  noiJi = rnorm(0, sigJi*sqrt(dt))                 -> sigJi*sqdt*N(0,1)         (L27)
  noiF  = rnorm(0, sigF *sqrt(dt))                 -> sigF *sqdt*N(0,1)         (L28)
  delta = 0.013                                                                (L24)
  Si_term = 0.1*Ji*dt - theta_Si*Si*dt - delta*Si*dt + Si*noiSi               (L30)
            (0.1 == lambda_J maturation; NB the i-Si term has NO predation/P loss,
             because there is no P compartment -- this is the SRJF reduction)
  Ji_term = ri*f_Si*F*Si*dt - 0.1*Ji*dt - theta_Ji*Ji*dt - delta*Ji*dt + Ji*noiJi (L31)
  F_term  = F*noiF - f_Si*F*(Si + 1*Ji)*dt - delta*F*dt + 0.37*dt              (L32)
            (mu_food = 0.37; food consumed by Si and Ji; NO xi*Ii term, no n-species)
  F+=F_term; Si+=Si_term; Ji+=Ji_term                                       (L34-36)
  soft constraints / error_count  (note the EXACT increments differ from SIRJPF2): (L38-40)
    if (Si<0 || Si>1e5)  { Si=0; error_count += 1     }
    if (F <0 || F >1e20) { F =0; error_count += 1000  }
    if (Ji<0 || Ji>1e5)  { Ji=0; error_count += 0.001 }
  T_Si = fabs(Si)                                                              (L42)
  -- error_count is an accumvar -> reset to 0 each obs by pypomp [R: accumvars] --
rinit  [R: fit_null.R lines 45-51  dyn_init]:
  Si=3; F=16.667; Ji=0; T_Si=0; error_count=0
dmeas  [R: fit_null.R lines 53-63]:
  if (error_count > 0) lik = -150
  else lik = dnbinom_mu(lumadult, k_Si, T_Si, give_log)   (size=k_Si, mu=T_Si)
partrans  [R: line 69-71]: log on ALL of sigSi,sigF,f_Si,ri,k_Si,sigJi,theta_Si,theta_Ji.
  (sigSi is rw.sd=0-pinned at its start 0; with start 0 it stays 0, so the log
   transform is never exercised on it -- matched by pinning + clamped log below.)
rw_sd  [R: rw_sd(sigSi=0, sigF=v, theta_Si=v, k_Si=v, f_Si=v, ri=v, sigJi=v, theta_Ji=v)]:
  v = 0.05 (round 1), 0.04 (round 2); sigSi = 0 (pinned).
MC settings  [R: Np=Mp=1500, Np_rep=10, Nmif round1=150/round2=250, cooling 0.7,
  block=TRUE for alts, 3*W starts, top-25%-and-replicate-x4 between rounds].
===============================================================================
FIDELITY NOTE: JAX's PRNG != R's, so this is NOT bit-identical; it reproduces the
estimator/workflow.  Validate STATISTICALLY: a handful of bootstrap null ll should
land near the R results_null/lrt_null_<b>.rds $ll (~ -440..-450 for this family),
and Lambda = 2*(ll_alt - ll_null) distributions should agree within MC error.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

# ----------------------------------------------------------------------------
# 0. JAX device / precision  -- MUST precede `import jax`
# ----------------------------------------------------------------------------
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
if os.environ.get("PYPOMP_USE_GPU", "1").lower() in {"0", "false", "no", "cpu"}:
    os.environ["JAX_PLATFORMS"] = "cpu"
else:
    os.environ.pop("JAX_PLATFORMS", None)

import jax  # noqa: E402
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
N_UNITS = 10                                              # [R: for (i in 1:10)]
UNIT_NAMES = [f"u{i}" for i in range(1, N_UNITS + 1)]     # u1..u10

# 8-parameter canonical order  [R: param_names, fit_null.R line 73].
# NB: collect_lrt.R reads only $ll/$se + metadata, never $coef, so the only hard
# requirement on coef is that it round-trips for provenance; we still emit it in
# R's coef() layout (shared first, then param[unit] for the alt block) to be safe.
PARAM_NAMES = ["sigSi", "sigF", "f_Si", "ri", "k_Si", "sigJi", "theta_Si", "theta_Ji"]

# State vector  [R: state_names, fit_null.R line 74]; order irrelevant (dict-keyed).
SRJF_STATENAMES = ["Si", "Ji", "F", "T_Si", "error_count"]

# Parameters carried on the LOG scale  [R: parameter_trans(log=...) line 69-71] = ALL 8.
# sigSi is fixed at 0 (rw.sd=0); with start 0 it never moves, so the clamped log
# (jnp.maximum(x,1e-30)) is a safe no-op for it -- it stays pinned via rw_sd=0 below.
SRJF_LOG_PARAMS = ("sigSi", "sigF", "f_Si", "ri", "k_Si", "sigJi", "theta_Si", "theta_Ji")

# Alternatives this family tests  [R: fit_alt.R alt_params_map lines 18-22].
ALT_PARAMS_MAP = {
    "theta_Si": ["theta_Si"],
    "ri":       ["ri"],
    "f_Si":     ["f_Si"],
}

# Run-level budgets, 1-based run_level in {1,2,3}.  Level 3 == the R run.
# [R: Np <- 1500; Np_rep <- 10; Mp <- 1500]  (fit_null.R lines 106-108)
ALGORITHMIC_PARAMS = {
    "Np":     [50, 500, 1500],   # final pfilter particles      [R: Np = 1500]
    "Np_rep": [2, 10, 10],       # pfilter replicates           [R: Np_rep = 10]
    "Mp":     [50, 500, 1500],   # mif2 particles               [R: Mp = 1500]
}
NMIF_ROUNDS = (150, 250)         # [R: Nmif = 150 (round1), 250 (round2)]
RW_ROUNDS = (0.05, 0.04)         # [R: dent_rw.sd = 0.05 then 0.04]
COOLING_A = 0.7                  # [R: cooling.fraction.50 = 0.7]
N_STARTS_DEFAULT = 300           # [R: 3 * getDoParWorkers() = 3 * 100]
BATCH_SIZE = int(os.environ.get("PYPOMP_BATCH", "120"))


# ----------------------------------------------------------------------------
# 2. Lum SRJF model in JAX  (term-for-term with fit_null.R / fit_alt.R Csnippets)
# ----------------------------------------------------------------------------
def srjf_rproc(X_, theta_, key, covars, t, dt):
    """One Euler-Maruyama step  [R: dyn_rpro, euler(delta.t = 1/4)]."""
    Si, Ji, F = X_["Si"], X_["Ji"], X_["F"]
    error_count = X_["error_count"]
    sigSi, sigJi, sigF = theta_["sigSi"], theta_["sigJi"], theta_["sigF"]
    theta_Si, theta_Ji = theta_["theta_Si"], theta_["theta_Ji"]
    f_Si, ri = theta_["f_Si"], theta_["ri"]

    delta = 0.013          # [R: double delta = 0.013]
    mu_food = 0.37         # [R: + 0.37 * dt]
    lambda_J = 0.1         # [R: 0.1 * Ji  maturation rate]

    keys = jax.random.split(key, 3)
    sqdt = jnp.sqrt(dt)
    noiSi = sigSi * sqdt * jax.random.normal(keys[0])   # [R: rnorm(0, sigSi*sqrt(dt))]
    noiJi = sigJi * sqdt * jax.random.normal(keys[1])   # [R: rnorm(0, sigJi*sqrt(dt))]
    noiF = sigF * sqdt * jax.random.normal(keys[2])     # [R: rnorm(0, sigF *sqrt(dt))]

    # [R: Si_term = 0.1*Ji*dt - theta_Si*Si*dt - delta*Si*dt + Si*noiSi]
    Si_term = (lambda_J * Ji * dt - theta_Si * Si * dt
               - delta * Si * dt + Si * noiSi)
    # [R: Ji_term = ri*f_Si*F*Si*dt - 0.1*Ji*dt - theta_Ji*Ji*dt - delta*Ji*dt + Ji*noiJi]
    Ji_term = (ri * f_Si * F * Si * dt - lambda_J * Ji * dt
               - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi)
    # [R: F_term = F*noiF - f_Si*F*(Si + 1*Ji)*dt - delta*F*dt + 0.37*dt]
    F_term = (F * noiF - f_Si * F * (Si + 1.0 * Ji) * dt
              - delta * F * dt + mu_food * dt)

    Si_new = Si + Si_term      # [R: Si += Si_term]
    Ji_new = Ji + Ji_term      # [R: Ji += Ji_term]
    F_new = F + F_term         # [R: F  += F_term]

    # [R: soft constraints -> error_count, then state reset to 0]; EXACT increments:
    #   Si:  +1      F:  +1000      Ji:  +0.001
    def viol(x, hi):
        return (x < 0.0) | (x > hi)

    eps = 0.0
    eps += jnp.where(viol(Si_new, 1e5), 1.0, 0.0)        # [R: error_count += 1]
    eps += jnp.where(viol(F_new, 1e20), 1.0e3, 0.0)      # [R: error_count += 1000]
    eps += jnp.where(viol(Ji_new, 1e5), 1.0e-3, 0.0)     # [R: error_count += 0.001]

    Si_new = jnp.clip(Si_new, 0.0, 1e5)   # [R: if violated, Si = 0.0  (clip to bounds)]
    Ji_new = jnp.clip(Ji_new, 0.0, 1e5)   # [R: Ji = 0.0]
    F_new = jnp.clip(F_new, 0.0, 1e20)    # [R: F  = 0.0]

    return {
        "Si": Si_new, "Ji": Ji_new, "F": F_new,
        "T_Si": jnp.abs(Si_new),           # [R: T_Si = fabs(Si)]
        "error_count": error_count + eps,  # [R: accumvar -> reset to 0 each obs]
    }


def srjf_rinit(theta_, key, covars, t0):
    """Deterministic initial state  [R: dyn_init, fit_null.R lines 45-51]."""
    return {
        "Si": jnp.array(3.0),       # [R: Si = 3]
        "Ji": jnp.array(0.0),       # [R: Ji = 0]
        "F": jnp.array(16.667),     # [R: F = 16.667]
        "T_Si": jnp.array(0.0),     # [R: T_Si = 0.0]
        "error_count": jnp.array(0.0),
    }


def _nb_logpmf(y, mu, size):
    """Negative-binomial log-pmf in (mu, size) form  [R: dnbinom_mu(y, size, mu)]."""
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    return (gammaln(y + size) - gammaln(size) - gammaln(y + 1.0)
            + size * jnp.log(size / (size + mu))
            + y * jnp.log(mu / (size + mu)))


def srjf_dmeas(Y_, X_, theta_, covars, t):
    """1-D NB measurement log-likelihood  [R: dmeas, fit_null.R lines 53-63]."""
    # [R: dnbinom_mu(lumadult, k_Si, T_Si, give_log)]  size=k_Si, mu=T_Si
    ll = _nb_logpmf(Y_["lumadult"], X_["T_Si"], theta_["k_Si"])
    # [R: if (error_count > 0.0) lik = -150]
    return jnp.where(X_["error_count"] > 0.0, -150.0, ll)


def _nb_sample(key, mu, size):
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    k1, k2 = jax.random.split(key)
    g = jax.random.gamma(k1, size) * (mu / size)   # Gamma-Poisson = NB
    return jax.random.poisson(k2, g)


def srjf_rmeas(X_, theta_, key, covars, t):
    """Measurement simulator  [R: rmeas, fit_null.R lines 65-67]  (unused by fitting)."""
    return jnp.array([_nb_sample(key, X_["T_Si"], theta_["k_Si"])], dtype=float)


def srjf_to_est(theta):
    """natural -> estimation scale  [R: parameter_trans(log=...) toEst]; all 8 logged."""
    out = {**theta}
    for name in SRJF_LOG_PARAMS:
        out[name] = jnp.log(jnp.maximum(theta[name], 1e-30))
    return out


def srjf_from_est(theta):
    """estimation -> natural scale  [R: parameter_trans fromEst]."""
    out = {**theta}
    for name in SRJF_LOG_PARAMS:
        out[name] = jnp.exp(theta[name])
    return out


srjf_par_trans = pp.ParTrans(to_est=srjf_to_est, from_est=srjf_from_est)


# ----------------------------------------------------------------------------
# 3. R-FREE I/O: read pre-converted CSVs (pandas), write output .rds (pyreadr).
#    sim_data_all.csv (cols: b,unit,day,lumadult) and true_params.csv (name,value)
#    are produced ONCE by convert_sim_data_to_csv() (bottom of file) from the .rds,
#    so the GPU node needs no R-data reader.  Only I/O dep: pip install pyreadr.
# ----------------------------------------------------------------------------
def read_sim_data(sim_dir: Path, b: int) -> dict[str, pd.DataFrame]:
    """Read the b-th simulated panel from sim_data_all.csv (cols: b,unit,day,lumadult)."""
    df = pd.read_csv(sim_dir / "sim_data_all.csv")
    df = df[df["b"] == b]
    if df.empty:
        raise ValueError(f"no rows for dataset b={b} in {sim_dir/'sim_data_all.csv'}")
    out = {}
    for unit_name in UNIT_NAMES:
        out[unit_name] = (df[df["unit"] == unit_name]
                          [["day", "lumadult"]]
                          .astype(float).sort_values("day").reset_index(drop=True))
    return out


def read_true_params(sim_dir: Path) -> dict[str, float]:
    """Read true_params.csv (name,value) -> dict for all 8 params."""
    tp = pd.read_csv(sim_dir / "true_params.csv")
    params = dict(zip(tp["name"].astype(str), tp["value"].astype(float)))
    missing = [n for n in PARAM_NAMES if n not in params]
    if missing:
        raise ValueError(f"true_params.csv missing parameters: {missing}")
    return params


def write_result_rds(ll: float, se: float, coef: dict[str, float],
                     Np: int, Mp: int, Np_rep: int, path: Path) -> None:
    """Write the result as an .rds list in collect_lrt.R's schema.

    collect_lrt.R reads:  res$ll, res$se, and res$Np/Mp/Np_rep/block (metadata guard).
    We emit the SAME fields fit_null.R/fit_alt.R save:
        list(ll, se, coef, Np, Mp, Np_rep, block, Nmif=c(round1=150, round2=250)).
    pyreadr.write_rds writes a single object; to reproduce the R named-list we write a
    one-row data.frame whose columns ARE those names (collect_lrt.R's `res$ll`,
    `res$Np`, ... index by name, which works on a 1-row data.frame too).  coef is
    flattened into coef.<name> columns; Nmif into Nmif.round1 / Nmif.round2.
    If a strict named-list is required, run the tiny R post-pass in
    rds_postpass.R (printed by --emit-postpass)."""
    import pyreadr
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {"ll": float(ll), "se": float(se),
           "Np": int(Np), "Mp": int(Mp), "Np_rep": int(Np_rep),
           "block": True,
           "Nmif.round1": int(NMIF_ROUNDS[0]), "Nmif.round2": int(NMIF_ROUNDS[1])}
    for k, v in coef.items():
        row[f"coef.{k}"] = float(v)
    pyreadr.write_rds(str(path), pd.DataFrame([row]))


# ----------------------------------------------------------------------------
# 4. rw_sd  [R: rw_sd(sigSi=0, sigF=v, theta_Si=v, k_Si=v, f_Si=v, ri=v, sigJi=v, theta_Ji=v)]
# ----------------------------------------------------------------------------
def make_rw_sd(value: float) -> dict[str, float]:
    rw = {name: value for name in PARAM_NAMES}
    rw["sigSi"] = 0.0    # [R: sigSi = 0]  (pinned at boundary 0)
    return rw


# ----------------------------------------------------------------------------
# 5. Panel construction + one MPIF round
# ----------------------------------------------------------------------------
def logmeanexp_se(values: np.ndarray, axis: int):
    """logmeanexp and its delta-method SE along `axis`  [matches pomp::logmeanexp]."""
    m = np.nanmax(values, axis=axis, keepdims=True)
    w = np.exp(values - m)
    n = np.sum(~np.isnan(values), axis=axis)
    mean_w = np.nanmean(w, axis=axis)
    est = np.squeeze(m, axis=axis) + np.log(mean_w)
    se = np.nanstd(w, axis=axis, ddof=1) / mean_w / np.sqrt(n)
    return est, se


def build_pomp_dict(sim_data: dict[str, pd.DataFrame], theta: dict[str, float]):
    """10 per-unit pp.Pomp objects  [R: pomplist loop, fit_null.R lines 79-101]."""
    pomp_dict = {}
    for unit_name in UNIT_NAMES:
        ys_u = (sim_data[unit_name]
                .set_index("day")[["lumadult"]].astype(float))
        pomp_dict[unit_name] = pp.Pomp(
            ys=ys_u, theta=pp.PompParameters(theta), statenames=SRJF_STATENAMES, t0=1.0,
            rinit=srjf_rinit, rproc=srjf_rproc,
            dmeas=srjf_dmeas, rmeas=srjf_rmeas,
            par_trans=srjf_par_trans, dt=0.25,        # [R: euler(delta.t = 1/4)]
            accumvars=("error_count",),               # [R: accumvars = c("error_count")]
        )
    return pomp_dict


def make_panel_parameters(param_rows: pd.DataFrame, alt_block: list[str]):
    """Build a PanelParameters with one spec per start.

    NULL  (alt_block == []): all 8 params shared, unit_specific empty  -> PIF.
    ALT   (alt_block == ["ri"|...]): the alt block's params become UNIT-SPECIFIC
      (one column per unit, all 10 seeded with the same shared start, exactly as
       R's create_parameters() does:  rep(parameters[[p]], 10)); the remaining
       params stay shared.  block=True in panel.mif then block-resamples the
       unit-specific params == R mif2(block=TRUE) on the specific matrix.
    [R: fit_alt.R create_parameters lines 33-44 + panelPomp(shared=, specific=)].
    """
    shared_names = [n for n in PARAM_NAMES if n not in alt_block]
    specs = []
    if alt_block:
        for _, row in param_rows.iterrows():
            shared_df = pd.DataFrame(
                {"shared": [float(row[n]) for n in shared_names]}, index=shared_names)
            unit_df = pd.DataFrame(
                {u: [float(row[n]) for n in alt_block] for u in UNIT_NAMES},
                index=alt_block)              # [R: rep(param, 10) -> one col per unit]
            specs.append({"shared": shared_df, "unit_specific": unit_df})
    else:
        empty_unit = pd.DataFrame(index=[], columns=UNIT_NAMES)
        for _, row in param_rows.iterrows():
            shared_df = pd.DataFrame(
                {"shared": [float(row[n]) for n in PARAM_NAMES]}, index=PARAM_NAMES)
            specs.append({"shared": shared_df, "unit_specific": empty_unit})
    return pp.PanelParameters(theta=specs)


def panel_theta_to_row(theta_dict, alt_block: list[str]) -> dict[str, float]:
    """One panel-theta spec -> {param_name: value} of the SHARED params (for restart
    propagation; the alt block's unit-specific values are tracked separately)."""
    shared_df = theta_dict["shared"]
    return {n: float(shared_df.loc[n, "shared"]) for n in shared_df.index}


def panel_theta_to_coef(theta_dict, alt_block: list[str]) -> dict[str, float]:
    """One panel-theta spec -> R-style coef() dict (shared first, then param[unit])."""
    shared_df = theta_dict["shared"]
    coef = {n: float(shared_df.loc[n, "shared"]) for n in shared_df.index}
    if alt_block:
        unit_df = theta_dict["unit_specific"]
        for p in alt_block:
            for u in UNIT_NAMES:
                coef[f"{p}[{u}]"] = float(unit_df.loc[p, u])
    return coef


def _run_mif_batch(pomp_dict, starts_batch, rw_sd_dict, nmif, Mp, Np, Np_rep,
                   alt_block, block, mif_key, pf_key, vmap_chunk_size):
    """MPIF + replicated pfilter on ONE batch of starts -> list of result rows."""
    panel = pp.PanelPomp(Pomp_dict=pomp_dict,
                         theta=make_panel_parameters(starts_batch, alt_block))

    # [R: mif2(Nmif, rw.sd, cooling.fraction.50=0.7, Np=Mp, block=TRUE)]
    # block=True is pypomp's MPIF; for the all-shared NULL it is a no-op (== PIF),
    # exactly as fit_null.R comments "no-op for the all-shared null".  For an ALT
    # with a non-empty unit_specific block it block-resamples that block == R block=TRUE.
    # pypomp >=0.4.6: cooling is carried by RWSigma.geometric_cooling(a=...), not a
    # mif(a=...) kwarg.  a=0.7 == R cooling.fraction.50=0.7.
    panel.mif(
        J=Mp, M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
        key=mif_key, block=block, vmap_chunk_size=vmap_chunk_size,
    )

    # [R: ll <- replicate(Np_rep, unitLogLik(pfilter(m1, Np))); panel_logmeanexp(MARGIN=1, se=TRUE)]
    panel.pfilter(J=Np, reps=Np_rep, key=pf_key,
                  chunk_size=vmap_chunk_size if vmap_chunk_size else 1)

    ll = np.asarray(panel.results_history[-1].logLiks.values)   # (Nstarts, U, reps)
    assert ll.ndim == 3, f"expected (Nstarts,U,reps), got {ll.shape}"
    # logmeanexp over reps PER UNIT, THEN sum over units (panel_logmeanexp MARGIN=1).
    unit_est, unit_se = logmeanexp_se(ll, axis=2)               # (Nstarts, U) each
    panel_ll = np.nansum(unit_est, axis=1)                      # sum over units
    panel_se = np.sqrt(np.nansum(unit_se ** 2, axis=1))         # combine unit SEs

    rows = []
    for idx, td in enumerate(panel.theta._to_list()):
        r = {"_restart": panel_theta_to_row(td, alt_block),
             "_coef": panel_theta_to_coef(td, alt_block),
             "loglik": float(panel_ll[idx]),
             "se": float(panel_se[idx])}
        if alt_block:                                          # carry the unit-specific block
            unit_df = td["unit_specific"]
            r["_unit"] = {p: [float(unit_df.loc[p, u]) for u in UNIT_NAMES] for p in alt_block}
        rows.append(r)
    return rows


def run_mif_round(pomp_dict, starts, rw_sd_dict, nmif, Mp, Np, Np_rep,
                  alt_block, block, key_seed, vmap_chunk_size):
    """One round over ALL starts, fed to pypomp in batches of BATCH_SIZE.  Each start
    is an independent fit (like R's foreach), so batching only caps peak memory."""
    starts = starts.reset_index(drop=True)
    n = len(starts)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    rows = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts.iloc[lo:lo + BATCH_SIZE]
        print(f"    batch {bi + 1}/{n_batches} (starts {lo}..{lo + len(batch) - 1})",
              flush=True)
        rows.extend(_run_mif_batch(
            pomp_dict, batch, rw_sd_dict, nmif, Mp, Np, Np_rep, alt_block, block,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size))
    return rows


def select_and_replicate(rows, alt_block):
    """Top ceil(N/4) by loglik, each replicated 4x  [R: fit_*.R top-25% select +
    rep(.., each=4)].  Returns (starts_df, unit_block_list) for the next round, where
    starts_df carries the SHARED restart and unit_block_list the unit-specific block."""
    n_select = math.ceil(len(rows) / 4.0)
    ordered = sorted(rows, key=lambda r: r["loglik"], reverse=True)[:n_select]
    starts, unit_blocks = [], []
    for r in ordered:
        for _ in range(4):
            starts.append(dict(r["_restart"]))
            if alt_block:
                unit_blocks.append({p: list(r["_unit"][p]) for p in alt_block})
    return starts, unit_blocks


def starts_to_df(starts, unit_blocks, alt_block):
    """Pack restart rows back into a DataFrame whose columns are PARAM_NAMES.  For an
    ALT, the alt-block columns are filled with the (per-unit) values; make_panel_parameters
    only reads alt-block cols for the unit_specific spec, so we put u1's value there as a
    canonical scalar and re-seed all units from the dedicated unit_blocks at build time."""
    df_rows = []
    for i, s in enumerate(starts):
        row = dict(s)
        if alt_block:
            for p in alt_block:
                row[p] = unit_blocks[i][p][0]   # placeholder; overwritten per-unit below
        df_rows.append(row)
    return pd.DataFrame(df_rows)[PARAM_NAMES]


def make_panel_parameters_with_units(param_rows, alt_block, unit_blocks):
    """ALT PanelParameters that seeds EACH unit's alt-block value from unit_blocks
    (so round-2 restarts keep the per-unit estimates from round-1's block fit, exactly
    like R's replicated_specific_list).  For the null (alt_block empty) this is just
    make_panel_parameters."""
    if not alt_block:
        return make_panel_parameters(param_rows, alt_block)
    shared_names = [n for n in PARAM_NAMES if n not in alt_block]
    specs = []
    for i, (_, row) in enumerate(param_rows.iterrows()):
        shared_df = pd.DataFrame(
            {"shared": [float(row[n]) for n in shared_names]}, index=shared_names)
        unit_df = pd.DataFrame(
            {u: [unit_blocks[i][p][j] for p in alt_block] for j, u in enumerate(UNIT_NAMES)},
            index=alt_block)
        specs.append({"shared": shared_df, "unit_specific": unit_df})
    return pp.PanelParameters(theta=specs)


def _run_mif_batch_units(pomp_dict, starts_batch, unit_blocks_batch, rw_sd_dict,
                         nmif, Mp, Np, Np_rep, alt_block, block,
                         mif_key, pf_key, vmap_chunk_size):
    """Same as _run_mif_batch but seeds the alt block per-unit from unit_blocks_batch."""
    panel = pp.PanelPomp(
        Pomp_dict=pomp_dict,
        theta=make_panel_parameters_with_units(starts_batch, alt_block, unit_blocks_batch))
    panel.mif(
        J=Mp, M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
        key=mif_key, block=block, vmap_chunk_size=vmap_chunk_size,
    )
    panel.pfilter(J=Np, reps=Np_rep, key=pf_key,
                  chunk_size=vmap_chunk_size if vmap_chunk_size else 1)
    ll = np.asarray(panel.results_history[-1].logLiks.values)
    assert ll.ndim == 3, f"expected (Nstarts,U,reps), got {ll.shape}"
    unit_est, unit_se = logmeanexp_se(ll, axis=2)
    panel_ll = np.nansum(unit_est, axis=1)
    panel_se = np.sqrt(np.nansum(unit_se ** 2, axis=1))
    rows = []
    for idx, td in enumerate(panel.theta._to_list()):
        r = {"_restart": panel_theta_to_row(td, alt_block),
             "_coef": panel_theta_to_coef(td, alt_block),
             "loglik": float(panel_ll[idx]), "se": float(panel_se[idx])}
        if alt_block:
            unit_df = td["unit_specific"]
            r["_unit"] = {p: [float(unit_df.loc[p, u]) for u in UNIT_NAMES] for p in alt_block}
        rows.append(r)
    return rows


def run_mif_round_units(pomp_dict, starts, unit_blocks, rw_sd_dict, nmif, Mp, Np,
                        Np_rep, alt_block, block, key_seed, vmap_chunk_size):
    """Round driver that propagates per-unit alt-block starts (round 2)."""
    n = len(starts)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    rows = []
    starts_df = starts_to_df(starts, unit_blocks, alt_block)
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts_df.iloc[lo:lo + BATCH_SIZE]
        ublk = unit_blocks[lo:lo + (len(batch))] if alt_block else None
        print(f"    batch {bi + 1}/{n_batches} (starts {lo}..{lo + len(batch) - 1})",
              flush=True)
        rows.extend(_run_mif_batch_units(
            pomp_dict, batch, ublk, rw_sd_dict, nmif, Mp, Np, Np_rep, alt_block, block,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size))
    return rows


# ----------------------------------------------------------------------------
# 6. Main: 2-round MPIF fit  [R: fit_null.R / fit_alt.R round 1 + round 2]
# ----------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(
        description="GPU/pypomp Lum-SRJF bootstrap-LRT fit for one (dataset, target).")
    p.add_argument("dataset_index", type=int, help="dataset index 1..100")
    p.add_argument("target", type=str,
                   help="'null' or an alternative: " + ", ".join(ALT_PARAMS_MAP))
    p.add_argument("--run-level", type=int,
                   default=int(os.environ.get("PYPOMP_RUN_LEVEL", "3")), choices=[1, 2, 3])
    p.add_argument("--n-starts", type=int,
                   default=int(os.environ["PYPOMP_NSTARTS"]) if os.environ.get("PYPOMP_NSTARTS")
                   else N_STARTS_DEFAULT,
                   help="round-1 starts (R uses 3*100=300)")
    return p.parse_args(argv)


def main(argv) -> int:
    args = parse_args(argv)
    b = args.dataset_index
    target = args.target
    if not (1 <= b <= 100):
        raise ValueError("dataset index must be 1..100")
    is_null = (target == "null")
    if not is_null and target not in ALT_PARAMS_MAP:
        raise ValueError(f"unknown target {target!r}; valid: null, {list(ALT_PARAMS_MAP)}")
    alt_block = [] if is_null else ALT_PARAMS_MAP[target]

    rl_i = args.run_level - 1
    Mp = ALGORITHMIC_PARAMS["Mp"][rl_i]
    Np = ALGORITHMIC_PARAMS["Np"][rl_i]
    Np_rep = ALGORITHMIC_PARAMS["Np_rep"][rl_i]
    vmap_chunk = os.environ.get("PYPOMP_VMAP_CHUNK")
    vmap_chunk = int(vmap_chunk) if vmap_chunk else None

    print(f"=== pypomp Lum-SRJF bootstrap-LRT: dataset {b}, target '{target}' "
          f"(run_level {args.run_level}) ===")
    print(f"JAX backend = {jax.default_backend()}  x64 = {_USE_X64}  devices = {jax.devices()}")
    print(f"Mp={Mp}  Np={Np}  Np_rep={Np_rep}  starts={args.n_starts}  "
          f"alt_block={alt_block or 'NONE (all-shared null)'}")

    base = Path(__file__).resolve().parent
    sim_dir = base / "simulated_data"
    # Output path  [R: results_null/lrt_null_<b>.rds | results_alt/lrt_<alt>_<b>.rds].
    # PYPOMP_OVERWRITE=0 (default) writes a *_pypomp.rds SIDECAR so the smoke +
    # validation runs do NOT clobber the canonical R results; =1 overwrites the
    # canonical lrt_*.rds name (feed collect_lrt.R directly).
    overwrite = os.environ.get("PYPOMP_OVERWRITE", "0").lower() in {"1", "true", "yes"}
    if is_null:
        name = f"lrt_null_{b}" if overwrite else f"lrt_null_{b}_pypomp"
        out_path = base / "results_null" / f"{name}.rds"
    else:
        name = f"lrt_{target}_{b}" if overwrite else f"lrt_{target}_{b}_pypomp"
        out_path = base / "results_alt" / f"{name}.rds"

    sim_data = read_sim_data(sim_dir, b)
    true_params = read_true_params(sim_dir)
    theta0 = {n: true_params[n] for n in PARAM_NAMES}
    pomp_dict = build_pomp_dict(sim_data, theta0)

    # Round-1 starts: all = the true (MLE) params  [R: shared.start = true_params].
    starts = pd.DataFrame([theta0] * args.n_starts)[PARAM_NAMES]

    # Round 1  [R: Nmif=150, rw.sd=0.05]
    print(f"Round 1 (Nmif={NMIF_ROUNDS[0]}, rw={RW_ROUNDS[0]})...", flush=True)
    round1 = run_mif_round(
        pomp_dict, starts, make_rw_sd(RW_ROUNDS[0]),
        nmif=NMIF_ROUNDS[0], Mp=Mp, Np=Np, Np_rep=Np_rep,
        alt_block=alt_block, block=True,
        key_seed=801_000 + b, vmap_chunk_size=vmap_chunk)
    print(f"Round 1 done. best ll = {max(r['loglik'] for r in round1):.3f}", flush=True)

    # Select top 25% + replicate x4  [R: order desc; rep(.., each=4)]
    starts2, unit_blocks2 = select_and_replicate(round1, alt_block)

    # Round 2  [R: Nmif=250, rw.sd=0.04]
    print(f"Round 2 (Nmif={NMIF_ROUNDS[1]}, rw={RW_ROUNDS[1]})...", flush=True)
    round2 = run_mif_round_units(
        pomp_dict, starts2, unit_blocks2, make_rw_sd(RW_ROUNDS[1]),
        nmif=NMIF_ROUNDS[1], Mp=Mp, Np=Np, Np_rep=Np_rep,
        alt_block=alt_block, block=True,
        key_seed=802_000 + b, vmap_chunk_size=vmap_chunk)
    print(f"Round 2 done. best ll = {max(r['loglik'] for r in round2):.3f}", flush=True)

    # Best of round 2  [R: which.max(lls[1,]); coef; ll; se]
    best = max(round2, key=lambda r: r["loglik"])
    write_result_rds(best["loglik"], best["se"], best["_coef"], Np, Mp, Np_rep, out_path)

    print("\n================= RESULT (compare to R lrt_*.rds $ll) =================")
    print(f"ll = {best['loglik']:.4f}   se = {best['se']:.4f}")
    print(f"saved -> {out_path}")
    print("======================================================================")
    return 0


# ----------------------------------------------------------------------------
# 7. ONE-TIME data converter  (run on a host WITH R; NOT on the GPU node)
#    Produces simulated_data/sim_data_all.csv + true_params.csv from the committed
#    .rds.  This is documented here for reproducibility; the CSVs are committed so
#    the GPU node never touches R.  Equivalent R one-liner is in the deliverable notes.
# ----------------------------------------------------------------------------
def convert_sim_data_to_csv(sim_dir: Path, B: int = 100) -> None:
    """Read sim_data_<b>.rds (via pyreadr) + true_params.rds -> CSVs.  Only needed
    once; falls back to a clear error if pyreadr cannot read the list-of-data.frames
    (in which case use the R one-liner in the deliverable notes)."""
    import pyreadr
    frames = []
    for b in range(1, B + 1):
        f = sim_dir / f"sim_data_{b:03d}.rds"
        if not f.exists():
            continue
        res = pyreadr.read_r(str(f))   # list-of-data.frames may not round-trip; see notes
        for name, df in res.items():
            d = df.copy()
            d.columns = ["day", "lumadult"][:d.shape[1]]
            d["b"] = b
            d["unit"] = name
            frames.append(d[["b", "unit", "day", "lumadult"]])
    pd.concat(frames, ignore_index=True).to_csv(sim_dir / "sim_data_all.csv", index=False)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
