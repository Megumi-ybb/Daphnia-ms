#!/usr/bin/env python3
"""
pypomp_bootstrap_lrt.py  (Single-species Dent / SRJF)
=====================================================
GPU / pypomp port of the parametric-bootstrap-LRT cluster computation for the
Single-species **Dent SRJF** family (Sn / Jn / F only; one observable
``dentadult``; 10 mesocosm units).  It replicates, ON GPU and WITHOUT R, the
two R workers used in the cluster bootstrap-LRT:

    fit_null.R   ->  python pypomp_bootstrap_lrt.py <b> null
    fit_alt.R    ->  python pypomp_bootstrap_lrt.py <b> <alt>     # alt in {theta_Sn, rn, f_Sn}

For dataset ``b`` it builds the 10-unit panelPomp from ``sim_data_<b>``, fits the
all-shared NULL (or the unit-specific ALTERNATIVE) by the SAME two-round
marginalized panel iterated filtering (MPIF) as the R scripts, pfilters the best
fit, and writes the result rds CONSUMED UNCHANGED by ``collect_lrt.R``:

    results_null/lrt_null_<b>.rds            (target = "null")
    results_alt/lrt_<alt>_<b>.rds            (target = an alternative name)

each a list(ll, se, coef, Np, Mp, Np_rep, block, Nmif=c(round1=150, round2=250)),
exactly the schema fit_null.R / fit_alt.R save (and collect_lrt.R reads).

Usage
-----
    python pypomp_bootstrap_lrt.py <b 1..100> <target>
      target = "null"                 -> all-shared MPIF fit  (fit_null.R)
      target = "theta_Sn"|"rn"|"f_Sn" -> unit-specific MPIF fit (fit_alt.R)

Environment toggles (defaults reproduce the R run):
    PYPOMP_USE_GPU = 1 (default) | 0/cpu       -> JAX device
    PYPOMP_X64     = 1 (default) | 0           -> float64 (R is double; keep 1)
    PYPOMP_RUN_LEVEL = 3 (default) | 1 | 2     -> particle/iteration budget
                                                  (3 == the R settings below)
    PYPOMP_BATCH      = 96 (default) | <int>   -> outer batch over the MPIF starts
                                                  (GPU-memory knob; no effect on results)
    PYPOMP_VMAP_CHUNK = (unset) | <int>        -> pypomp per-UNIT vmap chunk

REQUIRES the development pypomp (>=0.4.6, with PanelPomp / PanelParameters /
RWSigma / ParTrans) and a CUDA JAX.  Inputs are the committed CSVs
(simulated_data/sim_data_all.csv, simulated_data/true_params.csv) read with
pandas, so NO R and NO rdata reader are needed on the GPU node.  The only output
dependency is ``pyreadr`` (to write the .rds in collect_lrt.R's schema):
``pip install pyreadr``.

By default (PYPOMP_OVERWRITE=0) the script writes a SIDECAR
results_{null,alt}/lrt_*_<b>_pypomp.rds that does NOT touch the canonical R
result; set PYPOMP_OVERWRITE=1 to overwrite the canonical lrt_*_<b>.rds consumed
by collect_lrt.R.

================================  FIDELITY NOTE  ===============================
JAX's PRNG is NOT R's, so the particle draws and the mif perturbations are NOT
bit-identical to a run of fit_null.R / fit_alt.R.  What this reproduces EXACTLY
is the *estimator and workflow*:

  * the Dent-SRJF process / measurement model term-for-term with the R Csnippets
    in fit_null.R lines 21-67 (== fit_alt.R lines 51-97);
  * the all-shared NULL (empty unit_specific) and, for an alternative, the
    SINGLE-parameter unit-specific block (that one param becomes per-unit, all
    others shared) -- exactly fit_alt.R's create_parameters();
  * the identical two-round MPIF: Nmif 150 / 250, Np=Mp=1500, Np_rep=10,
    cooling.fraction.50=0.7 (geometric), block=TRUE, rw.sd 0.05 then 0.04 with
    sigSn pinned at 0; round-2 starts = top 25% of round 1, each replicated x4;
  * the final best-loglik selection and the panel_logmeanexp(MARGIN=1, se=TRUE)
    loglik+se (logmeanexp over reps PER UNIT, then summed over units).

So the two runs must be compared statistically (ll within Monte-Carlo error),
NOT byte-for-byte.  The R<->pypomp mapping is annotated inline as [R: <desc>].
===============================================================================

R-Csnippet -> JAX mapping for the Dent-SRJF rproc  [R: fit_null.R lines 21-43]
------------------------------------------------------------------------------
R:  noiSn = rnorm(0, sigSn*sqrt(dt));  noiJn = rnorm(0, sigJn*sqrt(dt));
    noiF  = rnorm(0, sigF*sqrt(dt));
JAX: noiSn = sigSn*sqrt(dt)*N(0,1); etc.  (3 normals; one key split into 3)

R:  Sn_term = 0.1*Jn*dt - theta_Sn*Sn*dt          - delta*Sn*dt + Sn*noiSn;
JAX: Sn_term = lambda_J*Jn*dt - theta_Sn*Sn*dt    - delta*Sn*dt + Sn*noiSn   (lambda_J=0.1)
R:  Jn_term = rn*f_Sn*F*Sn*dt - 0.1*Jn*dt - theta_Jn*Jn*dt - delta*Jn*dt + Jn*noiJn;
JAX: identical (lambda_J=0.1)
R:  F_term  = F*noiF - f_Sn*F*(Sn + 1*Jn)*dt - delta*F*dt + 0.37*dt;
JAX: F_term = -f_Sn*F*(Sn + Jn)*dt - delta*F*dt + mu_food*dt + F*noiF   (mu_food=0.37)
R:  F+=F_term; Sn+=Sn_term; Jn+=Jn_term;   (Euler step, delta.t=1/4 -> dt=0.25)
R soft constraints (each sets the state to EXACTLY 0.0, not a clip-to-hi):
    if (Sn<0 || Sn>1e5){ Sn=0; error_count+=1; }
    if (F <0 || F >1e20){ F =0; error_count+=1000; }
    if (Jn<0 || Jn>1e5){ Jn=0; error_count+=0.001; }
JAX: viol(x,hi) = (x<0)|(x>hi); on viol set x->0.0 (matches R) and add the eps.
R:  T_Sn = fabs(Sn);                          JAX: T_Sn = abs(Sn_new).
error_count is an accumvar -> pypomp resets it to 0 at each observation.

rinit  [R: fit_null.R lines 45-51]:  Sn=3, F=16.667, Jn=0, T_Sn=0, error_count=0.
  (NB the all_shared model uses Sn=3 here, not 2.333 as in the SIRJPF2 port.)
dmeas  [R: fit_null.R lines 53-63]:  if (error_count>0) lik=-150 else
  lik = dnbinom_mu(dentadult, size=k_Sn, mu=T_Sn, log).   (ONE NB term)
rmeas  [R: fit_null.R lines 65-67]:  dentadult = rnbinom_mu(k_Sn, T_Sn).
partrans [R: fit_null.R lines 69-72]: log over
  (sigSn, sigF, f_Sn, rn, k_Sn, sigJn, theta_Sn, theta_Jn).
  sigSn is pinned at 0 (rw.sd=0); log(0)=-inf, so we treat sigSn as a FIXED-ZERO
  param (carried on the natural scale, never perturbed) -- numerically identical
  to R, where rw.sd=0 keeps sigSn at its boundary value 0.

DROPPED vs the SIRJPF2 JAX model (why):
  * states  In, Ii, Si, Ji, P, T_In, T_Si, T_Ii  -- this single-species n-only
    model has no infected (I), no invasive species (i/Si/Ii/Ji), no parasite (P);
  * params  sigIn,sigIi,sigSi,sigJi,sigP, theta_In,theta_Ii,theta_Si,theta_Ji,
    theta_P, f_Si, ri, probn, probi, xi, k_In,k_Si,k_Ii  -- all gone (no I/P/i);
  * the P+=25 inoculation pulse and its P>1e20 soft constraint (no parasite);
  * 3 of the 4 NB observables (dentinf, lumadult, luminf) -- only dentadult here.
  Kept: Sn, Jn, F (+ T_Sn, error_count) and params
    f_Sn, sigF, sigSn, rn, k_Sn, sigJn, theta_Sn, theta_Jn.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time
from pathlib import Path

# ----------------------------------------------------------------------------
# 0. JAX device / precision  -- MUST precede `import jax`
# ----------------------------------------------------------------------------
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")

_USE_GPU = os.environ.get("PYPOMP_USE_GPU", "1").lower() not in {"0", "false", "no", "cpu"}
if _USE_GPU:
    os.environ.pop("JAX_PLATFORMS", None)        # let JAX pick the GPU if present
else:
    os.environ["JAX_PLATFORMS"] = "cpu"

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
N_UNITS = 10                                            # [R: for i in 1:10]
UNIT_NAMES = [f"u{i}" for i in range(1, N_UNITS + 1)]   # u1..u10

# Canonical 8-parameter order  [R: param_names  fit_null.R line 74]
PARAM_NAMES = [
    "f_Sn", "sigF", "sigSn", "rn", "k_Sn", "sigJn", "theta_Sn", "theta_Jn",
]

# State vector  [R: state_names  fit_null.R line 76]; order irrelevant (dict-keyed).
DENT_STATENAMES = ["Sn", "Jn", "error_count", "F", "T_Sn"]

# Parameters carried on the LOG scale  [R: parameter_trans(log=...) fit_null.R 69-72].
# sigSn is in R's log list BUT is fixed at 0 (rw.sd=0); log(0)=-inf, so it is held
# as a FIXED-ZERO natural-scale param here, exactly matching R's pinned boundary value.
DENT_LOG_PARAMS = ("f_Sn", "sigF", "rn", "k_Sn", "sigJn", "theta_Sn", "theta_Jn")
FIXED_ZERO_PARAMS = ("sigSn",)                          # [R: sigSn=0, rw.sd=0]

# Alternatives this family tests (each makes ONE param unit-specific)  [R: fit_alt.R 18-22]
ALT_PARAMS_MAP = {
    "theta_Sn": ["theta_Sn"],
    "rn":       ["rn"],
    "f_Sn":     ["f_Sn"],
}

# --- MC settings  [R: fit_null.R/fit_alt.R lines 108-110, 113-148, 188] -----
# run_level 3 == the R cluster settings; 1/2 are smoke/cheap budgets only.
ALGORITHMIC_PARAMS = {
    "Np": [50, 500, 1500],      # final pfilter particles   [R: Np = 1500]
    "Np_rep": [2, 5, 10],       # pfilter replicates        [R: Np_rep = 10]
    "Mp": [50, 500, 1500],      # mif2 particles            [R: Mp = 1500]
}
NMIF_ROUND1 = 150               # [R: Nmif = 150  (round 1)]
NMIF_ROUND2 = 250               # [R: Nmif = 250  (round 2)]
RW_ROUND1 = 0.05                # [R: dent_rw.sd = 0.05  (round 1)]
RW_ROUND2 = 0.04                # [R: dent_rw.sd = 0.04  (round 2)]
COOLING_A = 0.7                 # [R: cooling.fraction.50 = 0.7]  (pypomp a == R)
N_STARTS_DEFAULT = int(os.environ.get("PYPOMP_NSTARTS", "300"))  # [R: 3*getDoParWorkers()=3*100=300]
                                # Each start is an independent fit (R foreach), so the
                                # count only trades compute for robustness; matches the
                                # cluster's 300 by default (other families use 300 too).
                                # Lower it (e.g. --n-starts 96 or PYPOMP_NSTARTS=96) only
                                # for a faster, lower-fidelity check.

# Outer batch over MPIF starts (GPU-memory knob; no effect on results -- each start
# is an independent fit, exactly like R's foreach).
BATCH_SIZE = int(os.environ.get("PYPOMP_BATCH", "96"))


# ----------------------------------------------------------------------------
# 2. Dent-SRJF model in JAX  [R: fit_null.R Csnippets lines 21-72]
#    (term-for-term with the R Csnippets; see the module docstring mapping)
# ----------------------------------------------------------------------------
def dent_rproc(X_, theta_, key, covars, t, dt):
    """One Euler-Maruyama step  [R: dyn_rpro, euler(delta.t = 1/4)]."""
    Sn, Jn, F = X_["Sn"], X_["Jn"], X_["F"]
    error_count = X_["error_count"]

    sigSn, sigJn, sigF = theta_["sigSn"], theta_["sigJn"], theta_["sigF"]
    theta_Sn, theta_Jn = theta_["theta_Sn"], theta_["theta_Jn"]
    f_Sn, rn = theta_["f_Sn"], theta_["rn"]

    delta = 0.013          # [R: double delta = 0.013]
    mu_food = 0.37         # [R: + 0.37 * dt]
    lambda_J = 0.1         # [R: 0.1 * Jn  (maturation rate)]

    keys = jax.random.split(key, 3)
    sqdt = jnp.sqrt(dt)
    noiSn = sigSn * sqdt * jax.random.normal(keys[0])   # [R: rnorm(0, sigSn*sqrt(dt))]
    noiJn = sigJn * sqdt * jax.random.normal(keys[1])   # [R: rnorm(0, sigJn*sqrt(dt))]
    noiF = sigF * sqdt * jax.random.normal(keys[2])     # [R: rnorm(0, sigF *sqrt(dt))]

    # [R: Sn_term = 0.1*Jn*dt - theta_Sn*Sn*dt - delta*Sn*dt + Sn*noiSn]
    Sn_term = (lambda_J * Jn * dt - theta_Sn * Sn * dt
               - delta * Sn * dt + Sn * noiSn)
    # [R: Jn_term = rn*f_Sn*F*Sn*dt - 0.1*Jn*dt - theta_Jn*Jn*dt - delta*Jn*dt + Jn*noiJn]
    Jn_term = (rn * f_Sn * F * Sn * dt - lambda_J * Jn * dt
               - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn)
    # [R: F_term = F*noiF - f_Sn*F*(Sn + 1*Jn)*dt - delta*F*dt + 0.37*dt]
    F_term = (-f_Sn * F * (Sn + Jn) * dt - delta * F * dt
              + mu_food * dt + F * noiF)

    Sn_new = Sn + Sn_term
    Jn_new = Jn + Jn_term
    F_new = F + F_term

    # [R: soft constraints -> error_count; the VIOLATING state is set to 0.0]
    def viol(x, hi):
        return (x < 0.0) | (x > hi)

    v_Sn = viol(Sn_new, 1e5)
    v_F = viol(F_new, 1e20)
    v_Jn = viol(Jn_new, 1e5)

    eps = 0.0
    eps += jnp.where(v_Sn, 1.0, 0.0)        # [R: error_count += 1]
    eps += jnp.where(v_F, 1.0e3, 0.0)       # [R: error_count += 1000]
    eps += jnp.where(v_Jn, 1.0e-3, 0.0)     # [R: error_count += 0.001]

    # [R: on violation the state is reset to EXACTLY 0.0 (NOT clipped to the bound)]
    Sn_new = jnp.where(v_Sn, 0.0, Sn_new)
    F_new = jnp.where(v_F, 0.0, F_new)
    Jn_new = jnp.where(v_Jn, 0.0, Jn_new)

    return {
        "Sn": Sn_new, "Jn": Jn_new, "F": F_new,
        "T_Sn": jnp.abs(Sn_new),            # [R: T_Sn = fabs(Sn)]
        "error_count": error_count + eps,   # accumvar -> reset to 0 at each obs by pypomp
    }


def dent_rinit(theta_, key, covars, t0):
    """Deterministic initial state  [R: dyn_init  fit_null.R lines 45-51]."""
    return {
        "Sn": jnp.array(3.0),               # [R: Sn = 3]
        "F": jnp.array(16.667),             # [R: F  = 16.667]
        "Jn": jnp.array(0.0),               # [R: Jn = 0]
        "T_Sn": jnp.array(0.0),             # [R: T_Sn = 0.0]
        "error_count": jnp.array(0.0),      # [R: error_count = 0.0]
    }


def _nb_logpmf(y, mu, size):
    """NB log-pmf in (mu, size) form  [R: dnbinom_mu(y, size=k_Sn, mu=T_Sn, log)]."""
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    return (gammaln(y + size) - gammaln(size) - gammaln(y + 1.0)
            + size * jnp.log(size / (size + mu))
            + y * jnp.log(mu / (size + mu)))


def dent_dmeas(Y_, X_, theta_, covars, t):
    """1-D NB measurement log-likelihood  [R: dmeas  fit_null.R lines 53-63]."""
    ll = _nb_logpmf(Y_["dentadult"], X_["T_Sn"], theta_["k_Sn"])
    # [R: if (error_count > 0.0) lik = -150]
    return jnp.where(X_["error_count"] > 0.0, -150.0, ll)


def _nb_sample(key, mu, size):
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    k1, k2 = jax.random.split(key)
    return jax.random.poisson(k2, jax.random.gamma(k1, size) * (mu / size))


def dent_rmeas(X_, theta_, key, covars, t):
    """Measurement simulator  [R: rmeas  fit_null.R lines 65-67]  (unused by the fit)."""
    return jnp.array([_nb_sample(key, X_["T_Sn"], theta_["k_Sn"])], dtype=float)


def dent_to_est(theta):
    """natural -> estimation scale  [R: parameter_trans(log=...) toEst]."""
    out = {**theta}
    for name in DENT_LOG_PARAMS:
        out[name] = jnp.log(jnp.maximum(theta[name], 1e-30))
    for name in FIXED_ZERO_PARAMS:
        out[name] = theta[name]            # sigSn pinned at 0, not log-transformed
    return out


def dent_from_est(theta):
    """estimation -> natural scale  [R: parameter_trans fromEst]."""
    out = {**theta}
    for name in DENT_LOG_PARAMS:
        out[name] = jnp.exp(theta[name])
    for name in FIXED_ZERO_PARAMS:
        out[name] = theta[name]
    return out


dent_par_trans = pp.ParTrans(to_est=dent_to_est, from_est=dent_from_est)


# ----------------------------------------------------------------------------
# 3. R-FREE I/O: read the pre-converted CSVs with pandas, write the .rds output
#    in collect_lrt.R's schema with pyreadr (required; raises clearly if absent --
#    no silent CSV/JSON fallback).
#    The CSVs (simulated_data/sim_data_all.csv, true_params.csv) are produced ONCE
#    by the R converter below from the real .rds files (which are NOT modified).
# ----------------------------------------------------------------------------
def read_sim_data(sim_dir: Path, b: int) -> dict[str, pd.DataFrame]:
    """Read the b-th simulated panel from sim_data_all.csv (cols: b,unit,day,dentadult)."""
    df = pd.read_csv(sim_dir / "sim_data_all.csv")
    df = df[df["b"] == b]
    if df.empty:
        raise ValueError(f"no rows for dataset b={b} in {sim_dir/'sim_data_all.csv'} "
                         f"(run the converter in this script's docstring first)")
    out = {}
    for unit_name in UNIT_NAMES:
        out[unit_name] = (df[df["unit"] == unit_name][["day", "dentadult"]]
                          .astype(float).sort_values("day").reset_index(drop=True))
    return out


def read_true_params(sim_dir: Path) -> dict[str, float]:
    """Read true_params.csv (name,value)  [R: simulated_data/true_params.rds]."""
    tp = pd.read_csv(sim_dir / "true_params.csv")
    params = dict(zip(tp["name"].astype(str), tp["value"].astype(float)))
    missing = [n for n in PARAM_NAMES if n not in params]
    if missing:
        raise ValueError(f"true_params.csv missing parameters: {missing}")
    return params


def write_result_rds(coef: dict[str, float], ll: float, se: float, Np: int, Mp: int,
                     Np_rep: int, out_path: Path) -> None:
    """Write list(ll, se, coef, Np, Mp, Np_rep, block, Nmif) as an .rds for
    collect_lrt.R, via pyreadr (EXACTLY as the SIRJPF2 reference port does).

    pyreadr cannot write an arbitrary R named-list, so we write a one-row
    data.frame whose columns hold every scalar field plus the coef vector spread
    across "coef.<name>" columns and Nmif.round1 / Nmif.round2.  collect_lrt.R
    reads res$ll, res$se, res$Np, res$Mp, res$Np_rep, res$block via `$`, which on
    a data.frame selects the column (a length-1 vector) -> identical value, so NO
    change to collect_lrt.R is needed for those fields.
    """
    import pyreadr
    out_path.parent.mkdir(parents=True, exist_ok=True)
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
    pyreadr.write_rds(str(out_path), pd.DataFrame(rec))


# ----------------------------------------------------------------------------
# 4. rw_sd  [R: rw_sd(sigSn=0, sigF=v, theta_Sn=v, k_Sn=v, f_Sn=v, rn=v, sigJn=v, theta_Jn=v)]
# ----------------------------------------------------------------------------
def make_rw_sd(value: float) -> dict[str, float]:
    """All 8 params perturbed at `value` EXCEPT sigSn (pinned at 0)  [R: fit_*.R rw.sd]."""
    rw = {name: value for name in PARAM_NAMES}
    rw["sigSn"] = 0.0                           # [R: sigSn = 0]
    return rw


# ----------------------------------------------------------------------------
# 5. Panel construction + one MPIF round (batched over starts for GPU memory)
# ----------------------------------------------------------------------------
def logmeanexp_se(values: np.ndarray, axis: int):
    """logmeanexp and its Monte-Carlo SE along `axis`  [matches pomp::logmeanexp(se=TRUE)]."""
    m = np.nanmax(values, axis=axis, keepdims=True)
    w = np.exp(values - m)
    n = np.sum(~np.isnan(values), axis=axis)
    mean_w = np.nanmean(w, axis=axis)
    est = np.squeeze(m, axis=axis) + np.log(mean_w)
    se = np.nanstd(w, axis=axis, ddof=1) / mean_w / np.sqrt(n)   # delta-method SE
    return est, se


def build_pomp_dict(sim_data: dict[str, pd.DataFrame], theta: dict[str, float]):
    """10 per-unit pp.Pomp objects  [R: pomplist loop  fit_null.R lines 81-103]."""
    pomp_dict = {}
    for unit_name in UNIT_NAMES:
        ys_u = (sim_data[unit_name].set_index("day")[["dentadult"]].astype(float))
        pomp_dict[unit_name] = pp.Pomp(
            ys=ys_u, theta=pp.PompParameters(theta), statenames=DENT_STATENAMES, t0=1.0,
            rinit=dent_rinit, rproc=dent_rproc, dmeas=dent_dmeas,
            rmeas=dent_rmeas, par_trans=dent_par_trans, dt=0.25,
            accumvars=("error_count",),         # [R: accumvars = c("error_count")]
        )
    return pomp_dict


def make_panel_parameters(param_rows: pd.DataFrame, specific_names: list[str]):
    """
    Build a PanelParameters with one spec per MPIF start.

    NULL  (specific_names == []):  all 8 params SHARED, empty unit_specific
        [R: panelPomp(pomplist, shared = true_params)].
    ALT   (specific_names == [p]): every param shared EXCEPT p, which becomes a
        per-unit row (one column per unit, all 10 seeded at the same value)
        [R: create_parameters(): shared = params[-p]; specific_mat = rep(p, 10)].
        panel.mif(block=True) is MPIF and block-resamples the unit-specific row.
    """
    shared_names = [n for n in PARAM_NAMES if n not in specific_names]
    specs = []
    for _, row in param_rows.iterrows():
        shared_df = pd.DataFrame({"shared": [float(row[n]) for n in shared_names]},
                                 index=shared_names)
        if specific_names:
            unit_df = pd.DataFrame(
                {u: [float(row[p]) for p in specific_names] for u in UNIT_NAMES},
                index=specific_names,
            )
        else:
            unit_df = pd.DataFrame(index=[], columns=UNIT_NAMES)
        specs.append({"shared": shared_df, "unit_specific": unit_df})
    return pp.PanelParameters(theta=specs)


def panel_theta_to_coef(theta_dict: dict, specific_names: list[str]) -> dict[str, float]:
    """
    Reconstruct an R-style coef() vector from one panel-theta spec:
    shared params by bare name, plus `p[u1]`..`p[u10]` for each unit-specific p.
    [R: coef(mif) layout, e.g. k_Sn, sigF, ..., rn[u1], ..., rn[u10]].
    """
    shared_df = theta_dict["shared"]
    coef = {n: float(shared_df.loc[n, "shared"]) for n in shared_df.index}
    if specific_names:
        unit_df = theta_dict["unit_specific"]
        for p in specific_names:
            for u in UNIT_NAMES:
                coef[f"{p}[{u}]"] = float(unit_df.loc[p, u])
    return coef


def panel_theta_to_shared_row(theta_dict: dict) -> dict[str, float]:
    """Flatten a fitted panel-theta's SHARED params to a {name: value} row, for the
    round-1 -> round-2 warm-start carry.  The unit-specific block is NOT collapsed
    here; it is carried separately, per-unit, into round 2 (see select_and_replicate
    and make_panel_parameters_with_units) -- matching fit_alt.R's
    replicated_specific_list, which keeps the FULL per-unit specific matrix."""
    shared_df = theta_dict["shared"]
    return {n: float(shared_df.loc[n, "shared"]) for n in shared_df.index}


def _run_mif_batch(pomp_dict, starts_batch, specific_names, rw_sd_dict, nmif,
                   Mp, Np, Np_rep, mif_key, pf_key, vmap_chunk_size):
    """MPIF + replicated pfilter on ONE batch of starts -> list of (coef, start, ll, se)."""
    panel = pp.PanelPomp(Pomp_dict=pomp_dict,
                         theta=make_panel_parameters(starts_batch, specific_names))

    # [R: mif2(Nmif, rw.sd, cooling.fraction.50=0.7, Np=Mp, block=TRUE)]
    # block=True is pypomp's MPIF; with a NON-EMPTY unit_specific (the alt) it
    # block-resamples that param, matching R's block=TRUE alt fit.  For the all-
    # shared null (empty unit_specific) block is a no-op == R's null block=TRUE.
    # pypomp >=0.4.6: cooling is carried by RWSigma.geometric_cooling(a=...),
    # not a mif(a=...) kwarg.  a=0.7 == R cooling.fraction.50=0.7.
    panel.mif(
        J=Mp, M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
        key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size,
    )

    # [R: ll <- replicate(Np_rep, unitLogLik(pfilter(m1, Np=Np)));
    #     panel_logmeanexp(ll, MARGIN=1, se=TRUE)]
    panel.pfilter(J=Np, reps=Np_rep, key=pf_key,
                  chunk_size=vmap_chunk_size if vmap_chunk_size else 1)

    # logLiks shape (Nstarts, U, reps).  logmeanexp over reps PER UNIT, then sum
    # over units (panel_logmeanexp(MARGIN=1)); SEs combined in quadrature.
    ll_array = np.asarray(panel.results_history[-1].logLiks.values)
    assert ll_array.ndim == 3, f"expected (Nstarts,U,reps), got {ll_array.shape}"
    unit_est, unit_se = logmeanexp_se(ll_array, axis=2)        # (Nstarts, U)
    panel_ll = np.nansum(unit_est, axis=1)                     # (Nstarts,)
    panel_se = np.sqrt(np.nansum(unit_se ** 2, axis=1))        # (Nstarts,)

    rows = []
    for idx, td in enumerate(panel.theta._to_list()):
        r = {
            "coef": panel_theta_to_coef(td, specific_names),
            "shared": panel_theta_to_shared_row(td),   # for round-2 warm-start carry
            "loglik": float(panel_ll[idx]),
            "se": float(panel_se[idx]),
        }
        if specific_names:                              # carry the FULL per-unit block
            unit_df = td["unit_specific"]
            r["unit"] = {p: [float(unit_df.loc[p, u]) for u in UNIT_NAMES]
                         for p in specific_names}
        rows.append(r)
    return rows


def run_mif_round(pomp_dict, starts, specific_names, rw_sd_dict, nmif,
                  Mp, Np, Np_rep, key_seed, vmap_chunk_size):
    """One round over ALL starts, fed to pypomp in batches of BATCH_SIZE.
    Each start is an independent fit (exactly like R's foreach), so batching does
    not change results -- it only caps peak GPU memory."""
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
            pomp_dict, batch, specific_names, rw_sd_dict, nmif, Mp, Np, Np_rep,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size))
    return rows


def select_and_replicate(rows: list[dict], specific_names: list[str]):
    """Top ceil(N/4) by loglik, each replicated x4, carrying BOTH the shared warm
    start AND (for the alt) the FULL per-unit unit-specific block from round 1.
    [R: fit_alt.R lines 176-188: select top 25%; shared rep(each=4); and the
    per-unit specific matrices kept verbatim in replicated_specific_list (NOT
    collapsed).]  Returns (starts, unit_blocks) aligned row-for-row, where
    starts is a list of shared-param dicts and unit_blocks a list of
    {param: [per-unit values]} dicts (empty for the null)."""
    n_select = math.ceil(len(rows) / 4.0)
    ordered = sorted(rows, key=lambda r: r["loglik"], reverse=True)[:n_select]
    starts, unit_blocks = [], []
    for r in ordered:
        for _ in range(4):                          # replicate x4
            starts.append(dict(r["shared"]))
            if specific_names:
                unit_blocks.append({p: list(r["unit"][p]) for p in specific_names})
    return starts, unit_blocks


def starts_to_df(starts: list[dict], unit_blocks: list[dict],
                 specific_names: list[str]) -> pd.DataFrame:
    """Pack the shared restart rows into a PARAM_NAMES-keyed DataFrame.  For an alt,
    the alt-block columns need SOME value (make_panel_parameters_with_units re-seeds
    each unit from unit_blocks, so this is only a placeholder); use u1's value."""
    df_rows = []
    for i, s in enumerate(starts):
        row = dict(s)
        if specific_names:
            for p in specific_names:
                row[p] = unit_blocks[i][p][0]       # placeholder; overwritten per-unit
        df_rows.append(row)
    return pd.DataFrame(df_rows)[PARAM_NAMES]


def make_panel_parameters_with_units(param_rows: pd.DataFrame, specific_names: list[str],
                                     unit_blocks: list[dict]):
    """Like make_panel_parameters, but the alt's unit_specific block is seeded
    PER UNIT from the carried round-1 estimates (unit_blocks[i][p][j]) rather than
    collapsed to one value -- so round 2 resumes the per-unit estimates exactly as
    fit_alt.R's specific.start = replicated_specific_list does.  For the null
    (empty specific_names) this delegates to make_panel_parameters."""
    if not specific_names:
        return make_panel_parameters(param_rows, specific_names)
    shared_names = [n for n in PARAM_NAMES if n not in specific_names]
    specs = []
    for i, (_, row) in enumerate(param_rows.iterrows()):
        shared_df = pd.DataFrame({"shared": [float(row[n]) for n in shared_names]},
                                 index=shared_names)
        unit_df = pd.DataFrame(
            {u: [float(unit_blocks[i][p][j]) for p in specific_names]
             for j, u in enumerate(UNIT_NAMES)},
            index=specific_names,
        )
        specs.append({"shared": shared_df, "unit_specific": unit_df})
    return pp.PanelParameters(theta=specs)


def _run_mif_batch_units(pomp_dict, starts_batch, unit_blocks_batch, specific_names,
                         rw_sd_dict, nmif, Mp, Np, Np_rep, mif_key, pf_key, vmap_chunk_size):
    """Same as _run_mif_batch but seeds the alt block per-unit from unit_blocks_batch
    (round-2 warm start carrying the full per-unit estimates)."""
    panel = pp.PanelPomp(
        Pomp_dict=pomp_dict,
        theta=make_panel_parameters_with_units(starts_batch, specific_names, unit_blocks_batch))
    panel.mif(
        J=Mp, M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
        key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size,
    )
    panel.pfilter(J=Np, reps=Np_rep, key=pf_key,
                  chunk_size=vmap_chunk_size if vmap_chunk_size else 1)
    ll_array = np.asarray(panel.results_history[-1].logLiks.values)
    assert ll_array.ndim == 3, f"expected (Nstarts,U,reps), got {ll_array.shape}"
    unit_est, unit_se = logmeanexp_se(ll_array, axis=2)
    panel_ll = np.nansum(unit_est, axis=1)
    panel_se = np.sqrt(np.nansum(unit_se ** 2, axis=1))
    rows = []
    for idx, td in enumerate(panel.theta._to_list()):
        rows.append({
            "coef": panel_theta_to_coef(td, specific_names),
            "loglik": float(panel_ll[idx]),
            "se": float(panel_se[idx]),
        })
    return rows


def run_mif_round_units(pomp_dict, starts, unit_blocks, specific_names, rw_sd_dict,
                        nmif, Mp, Np, Np_rep, key_seed, vmap_chunk_size):
    """Round-2 driver that propagates the per-unit alt-block starts (and shared
    starts) in batches, like run_mif_round but with the per-unit warm start."""
    n = len(starts)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    starts_df = starts_to_df(starts, unit_blocks, specific_names)
    rows = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts_df.iloc[lo:lo + BATCH_SIZE]
        ublk = unit_blocks[lo:lo + len(batch)] if specific_names else None
        print(f"    batch {bi + 1}/{n_batches} (starts {lo}..{lo + len(batch) - 1})",
              flush=True)
        rows.extend(_run_mif_batch_units(
            pomp_dict, batch, ublk, specific_names, rw_sd_dict, nmif, Mp, Np, Np_rep,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size))
    return rows


# ----------------------------------------------------------------------------
# 6. Main: one (dataset b, target) fit  [R: fit_null.R / fit_alt.R]
# ----------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(
        description="GPU/pypomp bootstrap-LRT fit for Dent-SRJF (one dataset, one target).")
    p.add_argument("b", type=int, help="dataset index 1..100")
    p.add_argument("target", type=str,
                   help="'null' (all-shared) or an alternative: "
                        + ", ".join(ALT_PARAMS_MAP))
    p.add_argument("--run-level", type=int,
                   default=int(os.environ.get("PYPOMP_RUN_LEVEL", "3")),
                   choices=[1, 2, 3], help="particle/iteration budget (default 3 == R)")
    p.add_argument("--n-starts", type=int, default=N_STARTS_DEFAULT,
                   help="number of MPIF starts (R cluster used 3*100=300)")
    return p.parse_args(argv)


def main(argv) -> int:
    args = parse_args(argv)
    b, target, rl = args.b, args.target, args.run_level
    if not (1 <= b <= 100):
        raise ValueError("dataset index b must be 1..100")
    is_null = (target == "null")
    if not is_null and target not in ALT_PARAMS_MAP:
        raise ValueError(f"unknown target {target!r}; valid: null, {list(ALT_PARAMS_MAP)}")
    specific_names = [] if is_null else ALT_PARAMS_MAP[target]

    rl_i = rl - 1
    Mp = ALGORITHMIC_PARAMS["Mp"][rl_i]
    Np = ALGORITHMIC_PARAMS["Np"][rl_i]
    Np_rep = ALGORITHMIC_PARAMS["Np_rep"][rl_i]
    vmap_chunk = os.environ.get("PYPOMP_VMAP_CHUNK")
    vmap_chunk = int(vmap_chunk) if vmap_chunk else None

    print(f"=== pypomp Dent-SRJF bootstrap-LRT: dataset {b}, target '{target}' "
          f"(run_level {rl}) ===")
    print(f"JAX backend = {jax.default_backend()}  x64 = {_USE_X64}  "
          f"devices = {jax.devices()}")
    print(f"Mp={Mp} Np={Np} Np_rep={Np_rep} starts={args.n_starts} batch={BATCH_SIZE}  "
          f"unit_specific={specific_names or '(none, all-shared null)'}")

    base = Path(__file__).resolve().parent
    sim_dir = base / "simulated_data"
    # Output path: PYPOMP_OVERWRITE=0 (default) writes a SIDECAR (*_pypomp.rds) that
    # does NOT touch the canonical R result; =1 overwrites the canonical lrt_*.rds
    # consumed by collect_lrt.R  [matches the SIRJPF2 reference port].
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

    # Round-1 starts: all = the true params  [R: shared.start = true_params; all
    # foreach workers start from the same MLE warm start].
    starts = pd.DataFrame([theta0] * args.n_starts)[PARAM_NAMES]

    # Seed offsets mirror R: null uses 801*1000+b; alt uses 801*1000+b+100000.
    seed_base = 801_000 + b + (0 if is_null else 100_000)

    print("Round 1 starting...", flush=True)
    _t0 = time.time()
    round1 = run_mif_round(pomp_dict, starts, specific_names, make_rw_sd(RW_ROUND1),
                           nmif=NMIF_ROUND1, Mp=Mp, Np=Np, Np_rep=Np_rep,
                           key_seed=seed_base + 1000, vmap_chunk_size=vmap_chunk)
    print(f"Round 1 done. best ll = {max(r['loglik'] for r in round1):.3f} "
          f"({(time.time() - _t0) / 60:.1f} min)", flush=True)

    # Select top 25% + replicate x4, CARRYING the full per-unit alt block (not
    # collapsed)  [R: fit_alt.R replicated_specific_list].
    starts2, unit_blocks2 = select_and_replicate(round1, specific_names)

    print("Round 2 starting...", flush=True)
    _t1 = time.time()
    round2 = run_mif_round_units(pomp_dict, starts2, unit_blocks2, specific_names,
                                 make_rw_sd(RW_ROUND2),
                                 nmif=NMIF_ROUND2, Mp=Mp, Np=Np, Np_rep=Np_rep,
                                 key_seed=seed_base + 2000, vmap_chunk_size=vmap_chunk)
    print(f"Round 2 done. ({(time.time() - _t1) / 60:.1f} min)  "
          f"[job total {(time.time() - _t0) / 60:.1f} min]", flush=True)

    # [R: best <- which.max(lls[1,]); result = list(ll, se, coef, ...)]
    best = max(round2, key=lambda r: r["loglik"])
    write_result_rds(best["coef"], best["loglik"], best["se"], Np, Mp, Np_rep, out_path)

    coef = best["coef"]
    print("\n================= RESULT =================")
    print(f"ll = {best['loglik']:.4f}   se = {best['se']:.4f}")
    print(f"coef entries = {len(coef)}  (shared {len([k for k in coef if '[' not in k])} "
          f"+ unit-specific {len([k for k in coef if '[' in k])})")
    print(f"saved -> {out_path}")
    print("==========================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
