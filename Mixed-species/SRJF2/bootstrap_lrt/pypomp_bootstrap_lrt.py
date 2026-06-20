#!/usr/bin/env python3
"""
pypomp_bootstrap_lrt.py  (SRJF2 / Mixed-species)
================================================
GPU / pypomp port of the SRJF2 parametric-bootstrap-LRT cluster fits
(fit_null.R + fit_alt.R).  It reproduces, for ONE bootstrap dataset b and ONE
target (the all-shared NULL, or one block-unit-specific ALTERNATIVE), the SAME
two-round marginalized panel iterated filtering (MPIF) fit + replicated pfilter
log-likelihood, and writes the SAME .rds schema that collect_lrt.R consumes
unchanged:

    results_null/lrt_null_<b>.rds                     (target = "null")
    results_alt/lrt_<alt>_<b>.rds                     (target = an alternative)

with fields list(ll, se, coef, Np, Mp, Np_rep, block, Nmif=c(round1=150, round2=250)).

Usage
-----
    python pypomp_bootstrap_lrt.py <b 1..100> <target>

    target = "null"                  -> all-shared fit (empty unit_specific)
           | "f_Si_f_Sn"             -> (f_Si, f_Sn) made UNIT-SPECIFIC, rest shared
           | "rn_ri"                 -> (rn, ri)    made UNIT-SPECIFIC, rest shared
           | "theta_Sn_theta_Si"     -> (theta_Sn, theta_Si) UNIT-SPECIFIC, rest shared

    examples:
        python pypomp_bootstrap_lrt.py 7 null
        python pypomp_bootstrap_lrt.py 7 rn_ri
        python pypomp_bootstrap_lrt.py 7 f_Si_f_Sn
        python pypomp_bootstrap_lrt.py 7 theta_Sn_theta_Si

Environment toggles (defaults reproduce the R run):
    PYPOMP_USE_GPU = 1 (default) | 0/cpu
    PYPOMP_X64     = 1 (default)                 # R/pomp likelihoods are double precision
    PYPOMP_RUN_LEVEL = 3 (default) | 1 | 2       # particle/iteration budget (3 == R settings)
    PYPOMP_NSTARTS  = (auto) | <int>             # number of search starts (R: 3*100=300)
    PYPOMP_BATCH    = 100 (default) | <int>      # outer batch over starts (GPU memory)
    PYPOMP_VMAP_CHUNK = (unset) | <int>          # pypomp per-unit vmap chunk
    PYPOMP_OVERWRITE = 0 (default) | 1           # 0 -> write *_pypomp.rds SIDECAR (do NOT
                                                 #      touch the canonical R lrt_*.rds);
                                                 # 1 -> overwrite the canonical lrt_*.rds

REQUIRES the development pypomp (>=0.4.6, with PanelPomp/PanelParameters/RWSigma/ParTrans)
and a CUDA JAX.  R-FREE on the GPU node: the inputs are the committed CSVs
(simulated_data/sim_data_all.csv, simulated_data/true_params.csv), read with pandas;
the only output-I/O dependency is `pip install pyreadr` (to write the .rds).

================================  FIDELITY NOTE  ===============================
JAX's PRNG is NOT R's, so this is NOT bit-identical to fit_null.R / fit_alt.R; it
reproduces the same ESTIMATOR + WORKFLOW (model, MPIF structure, particle/iteration
counts, cooling, rw sds, top-25%-and-replicate selection, block=TRUE MPIF for the
alts).  Compare the per-b ll STATISTICALLY (agree within Monte-Carlo error), not
byte-for-byte.  The R<->pypomp mapping is annotated inline as  [R: <line/desc>].

------------------------------  MODEL DERIVATION  -----------------------------
The SRJF JAX model below is the SIRJPF2 model (pypomp_all_shared.py) REDUCED to
the SRJF Csnippet in fit_null.R / fit_alt.R (lines 21-90).  Dropped vs SIRJPF2:

  * COMPARTMENTS dropped: In, Ii (infected), P (parasite).  SRJF keeps only
    Sn, Jn (native S/juvenile), Si, Ji (invasive S/juvenile), F (food).
  * PARAMS dropped: xi, probn, probi, theta_In, theta_Ii, theta_P, sigP,
    sigIn, sigIi, k_In, k_Ii.  (SRJF2 has 15 params, SIRJPF2 has 26.)
  * OBSERVABLES dropped: dentinf, luminf.  SRJF measures only the two adult
    counts dentadult (~T_Sn) and lumadult (~T_Si).
  * The P+=25 inoculation pulse at day 4 is dropped (no parasite).
  * No xi-weighting in the F equation: F loses food to (Sn + Jn) and (Si + Ji)
    only (the SIRJPF2 (Sn + xi*In + Jn) collapses to (Sn + Jn) when In is gone).
  Kept identical: dt=0.25 euler, delta=0.013, lambda_J=0.1 maturation, mu_food=0.37,
  the NB(mu,size) measurement, the soft-constraint error_count + clip-to-[0,hi],
  T_Sn=fabs(Sn) / T_Si=fabs(Si), and the log parameter transform.
  See the per-term [R: ...] annotations in srjf_rproc / srjf_rinit / srjf_dmeas.
===============================================================================
"""

from __future__ import annotations

import argparse
import math
import os
import sys
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

# float64: R/pomp likelihoods are double precision; keep x64 ON for fidelity.
_USE_X64 = os.environ.get("PYPOMP_X64", "1").lower() not in {"0", "false", "no"}
jax.config.update("jax_enable_x64", _USE_X64)

import jax.numpy as jnp  # noqa: E402
from jax.scipy.special import gammaln  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pypomp as pp  # noqa: E402


# ----------------------------------------------------------------------------
# 1. Constants  [R: fit_null.R / fit_alt.R lines 88-118]
# ----------------------------------------------------------------------------
N_UNITS = 9                                          # [R: for (i in 1:9)]
UNIT_NAMES = [f"u{i}" for i in range(1, N_UNITS + 1)]

# Run-level budgets, indexed by run_level in {1,2,3} (1-based, matching R).  The
# R cluster fit uses Np=Mp=1500, Np_rep=10 (== run_level 3).  Levels 1/2 are
# CPU-smoke budgets ONLY (not used for the real GPU run).
# [R: fit_null.R lines 126-128  Np=1500; Np_rep=10; Mp=1500]
ALGORITHMIC_PARAMS = {
    "Np": [20, 500, 1500],       # final pfilter particles      [R: Np = 1500]
    "Np_rep": [2, 10, 10],       # pfilter replicates           [R: Np_rep = 10]
    "Mp": [20, 500, 1500],       # mif2 particles               [R: Mp = 1500]
}
NMIF_ROUND1 = 150               # [R: fit_*.R  Nmif = 150  (round 1)]
NMIF_ROUND2 = 250               # [R: fit_*.R  Nmif = 250  (round 2)]
RW_ROUND1 = 0.05                # [R: dent_rw.sd = 0.05  (round 1)]
RW_ROUND2 = 0.04                # [R: dent_rw.sd = 0.04  (round 2)]
COOLING_A = 0.7                 # [R: cooling.fraction.50 = 0.7]  (pypomp a == R cooling.fraction.50)

# Number of search starts.  R uses 3 * getDoParWorkers() = 3 * 100 = 300 starts.
# [R: foreach(i = 1:(3 * getDoParWorkers())) with registerDoParallel(cores = 100)]
N_STARTS_DEFAULT = int(os.environ.get("PYPOMP_NSTARTS", "300"))

# Outer batch size over STARTS (the GPU-memory / pass-count knob).  Each start is
# an independent fit (exactly like R's foreach), so batching does NOT change
# results -- only peak memory / number of sequential passes.
BATCH_SIZE = int(os.environ.get("PYPOMP_BATCH", "100"))

# Canonical 15-parameter order  [R: param_names  fit_null.R lines 88-90]
PARAM_NAMES = [
    "sigSn", "sigSi", "sigF", "f_Sn", "f_Si",
    "rn", "ri", "k_Sn", "k_Si", "sigJi", "sigJn",
    "theta_Sn", "theta_Si", "theta_Jn", "theta_Ji",
]

# State vector  [R: state_names  fit_null.R line 91];  order irrelevant (dict-keyed).
SRJF_STATENAMES = ["Sn", "Si", "Jn", "Ji", "error_count", "F", "T_Sn", "T_Si"]

# Parameters carried on the LOG scale  [R: parameter_trans(log=...)  fit_null.R lines 82-86].
# sigSn, sigSi are fixed at 0 and are NOT log-transformed here (log(0)=-inf); they
# are pinned by rw_sd = 0, exactly as in R where rw.sd = 0 keeps them at 0.
SRJF_LOG_PARAMS = (
    "sigF", "f_Sn", "f_Si", "rn", "ri", "k_Sn", "k_Si",
    "sigJi", "sigJn", "theta_Sn", "theta_Si", "theta_Jn", "theta_Ji",
)
FIXED_ZERO_PARAMS = ("sigSn", "sigSi")            # [R: rw_sd(sigSn=0, sigSi=0)]

# The three alternatives this family tests -> the block of params made UNIT-SPECIFIC.
# [R: fit_alt.R alt_params_map lines 18-22]
ALT_PARAMS_MAP = {
    "f_Si_f_Sn":         ["f_Si", "f_Sn"],
    "rn_ri":             ["rn", "ri"],
    "theta_Sn_theta_Si": ["theta_Sn", "theta_Si"],
}


# ----------------------------------------------------------------------------
# 2. SRJF model in JAX  [R: Csnippets  fit_null.R / fit_alt.R lines 21-90]
#    Derived TERM-BY-TERM from the SRJF Csnippet (see MODEL DERIVATION header).
# ----------------------------------------------------------------------------
def srjf_rproc(X_, theta_, key, covars, t, dt):
    """One Euler-Maruyama step  [R: dyn_rpro, euler(delta.t = 1/4)  lines 21-52]."""
    Sn, Jn = X_["Sn"], X_["Jn"]
    Si, Ji = X_["Si"], X_["Ji"]
    F = X_["F"]
    error_count = X_["error_count"]

    sigSn, sigSi = theta_["sigSn"], theta_["sigSi"]
    sigJn, sigJi = theta_["sigJn"], theta_["sigJi"]
    sigF = theta_["sigF"]
    theta_Sn, theta_Si = theta_["theta_Sn"], theta_["theta_Si"]
    theta_Jn, theta_Ji = theta_["theta_Jn"], theta_["theta_Ji"]
    f_Sn, f_Si = theta_["f_Sn"], theta_["f_Si"]
    rn, ri = theta_["rn"], theta_["ri"]

    delta = 0.013          # [R: double delta = 0.013]
    mu_food = 0.37         # [R: + 0.37 * dt]
    lambda_J = 0.1         # [R: 0.1 * Jn ... maturation rate]

    # [R: noiSn = rnorm(0, sigSn*sqrt(dt)); etc.]  5 process-noise draws (no In/Ii/P noise).
    keys = jax.random.split(key, 5)
    sqdt = jnp.sqrt(dt)
    noiSn = sigSn * sqdt * jax.random.normal(keys[0])
    noiSi = sigSi * sqdt * jax.random.normal(keys[1])
    noiJn = sigJn * sqdt * jax.random.normal(keys[2])
    noiJi = sigJi * sqdt * jax.random.normal(keys[3])
    noiF = sigF * sqdt * jax.random.normal(keys[4])

    # [R: Sn_term = 0.1*Jn*dt - theta_Sn*Sn*dt - delta*Sn*dt + Sn*noiSn]
    Sn_term = (lambda_J * Jn * dt - theta_Sn * Sn * dt
               - delta * Sn * dt + Sn * noiSn)
    # [R: Jn_term = rn*f_Sn*F*Sn*dt - 0.1*Jn*dt - theta_Jn*Jn*dt - delta*Jn*dt + Jn*noiJn]
    Jn_term = (rn * f_Sn * F * Sn * dt - lambda_J * Jn * dt
               - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn)
    # [R: Si_term = 0.1*Ji*dt - theta_Si*Si*dt - delta*Si*dt + Si*noiSi]
    Si_term = (lambda_J * Ji * dt - theta_Si * Si * dt
               - delta * Si * dt + Si * noiSi)
    # [R: Ji_term = ri*f_Si*F*Si*dt - 0.1*Ji*dt - theta_Ji*Ji*dt - delta*Ji*dt + Ji*noiJi]
    Ji_term = (ri * f_Si * F * Si * dt - lambda_J * Ji * dt
               - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi)
    # [R: F_term = F*noiF - f_Sn*F*(Sn + 1*Jn)*dt - f_Si*F*(Si + 1*Ji)*dt - delta*F*dt + 0.37*dt]
    F_term = (F * noiF
              - f_Sn * F * (Sn + Jn) * dt
              - f_Si * F * (Si + Ji) * dt
              - delta * F * dt + mu_food * dt)

    # [R: F += F_term; Sn += Sn_term; Si += Si_term; Ji += Ji_term; Jn += Jn_term]
    Sn_new = Sn + Sn_term
    Jn_new = Jn + Jn_term
    Si_new = Si + Si_term
    Ji_new = Ji + Ji_term
    F_new = F + F_term

    # [R: soft constraints -> error_count, then offending state -> 0]
    #   Sn:  err += 1            (line 44)
    #   Si:  err += 1000000      (line 45)
    #   F :  err += 1000         (line 46)
    #   Jn:  err += 0.001        (line 47)
    #   Ji:  err += 0.000000001  (line 48)
    def viol(x, hi):
        return (x < 0.0) | (x > hi)

    eps = 0.0
    eps += jnp.where(viol(Sn_new, 1e5), 1.0, 0.0)
    eps += jnp.where(viol(Si_new, 1e5), 1.0e6, 0.0)
    eps += jnp.where(viol(F_new, 1e20), 1.0e3, 0.0)
    eps += jnp.where(viol(Jn_new, 1e5), 1.0e-3, 0.0)
    eps += jnp.where(viol(Ji_new, 1e5), 1.0e-9, 0.0)

    # [R: in C the offending state is set to 0; pypomp clips to [0, hi] which is
    #  equivalent for the kept (non-violating) particles and bounded for the rest.
    #  We clip-to-0 the lower bound (matches "= 0.0") and to hi the upper bound.]
    Sn_new = jnp.clip(Sn_new, 0.0, 1e5)
    Jn_new = jnp.clip(Jn_new, 0.0, 1e5)
    Si_new = jnp.clip(Si_new, 0.0, 1e5)
    Ji_new = jnp.clip(Ji_new, 0.0, 1e5)
    F_new = jnp.clip(F_new, 0.0, 1e20)

    return {
        "Sn": Sn_new, "Jn": Jn_new, "Si": Si_new, "Ji": Ji_new, "F": F_new,
        # [R: T_Sn = fabs(Sn); T_Si = fabs(Si)]
        "T_Sn": jnp.abs(Sn_new), "T_Si": jnp.abs(Si_new),
        # [R: error_count is an accumvar -> reset to 0 at each obs by pypomp]
        "error_count": error_count + eps,
    }


def srjf_rinit(theta_, key, covars, t0):
    """Deterministic initial state  [R: dyn_init  fit_null.R lines 54-63]."""
    return {
        "Sn": jnp.array(2.333), "Si": jnp.array(0.667),     # [R: Sn=2.333; Si=0.667]
        "Jn": jnp.array(0.0), "Ji": jnp.array(0.0),          # [R: Jn=0; Ji=0]
        "F": jnp.array(16.667),                              # [R: F=16.667]
        "T_Sn": jnp.array(0.0), "T_Si": jnp.array(0.0),      # [R: T_Sn=0; T_Si=0]
        "error_count": jnp.array(0.0),                       # [R: error_count=0]
    }


def _srjf_nb_logpmf(y, mu, size):
    """Negative-binomial log-pmf in (mu, size) form  [R: dnbinom_mu(y, size, mu)]."""
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    return (gammaln(y + size) - gammaln(size) - gammaln(y + 1.0)
            + size * jnp.log(size / (size + mu))
            + y * jnp.log(mu / (size + mu)))


def srjf_dmeas(Y_, X_, theta_, covars, t):
    """2-D NB measurement log-likelihood  [R: dmeas  fit_null.R lines 65-75].
    [R: lik = dnbinom_mu(lumadult, k_Si, T_Si, log) + dnbinom_mu(dentadult, k_Sn, T_Sn, log)]
    (size = k_*, mu = T_*).  Note: only the two ADULT counts are measured."""
    error_count = X_["error_count"]
    ll = (_srjf_nb_logpmf(Y_["dentadult"], X_["T_Sn"], theta_["k_Sn"])
          + _srjf_nb_logpmf(Y_["lumadult"], X_["T_Si"], theta_["k_Si"]))
    # [R: if (error_count > 0.0) lik = -150]
    return jnp.where(error_count > 0.0, -150.0, ll)


def _srjf_nb_sample(key, mu, size):
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    k1, k2 = jax.random.split(key)
    g = jax.random.gamma(k1, size) * (mu / size)   # Gamma-Poisson = NB
    return jax.random.poisson(k2, g)


def srjf_rmeas(X_, theta_, key, covars, t):
    """Measurement simulator  [R: rmeas  fit_null.R lines 77-80]  (unused by fitting)."""
    keys = jax.random.split(key, 2)
    return jnp.array([
        _srjf_nb_sample(keys[0], X_["T_Sn"], theta_["k_Sn"]),   # dentadult
        _srjf_nb_sample(keys[1], X_["T_Si"], theta_["k_Si"]),   # lumadult
    ], dtype=float)


def srjf_to_est(theta):
    """natural -> estimation scale  [R: parameter_trans(log=...) toEst]."""
    out = {**theta}
    for name in SRJF_LOG_PARAMS:
        out[name] = jnp.log(jnp.maximum(theta[name], 1e-30))
    for name in FIXED_ZERO_PARAMS:
        out[name] = theta[name]        # pinned at 0, not log-transformed
    return out


def srjf_from_est(theta):
    """estimation -> natural scale  [R: parameter_trans fromEst]."""
    out = {**theta}
    for name in SRJF_LOG_PARAMS:
        out[name] = jnp.exp(theta[name])
    for name in FIXED_ZERO_PARAMS:
        out[name] = theta[name]
    return out


srjf_par_trans = pp.ParTrans(to_est=srjf_to_est, from_est=srjf_from_est)


# ----------------------------------------------------------------------------
# 3. R-FREE I/O: read the pre-converted CSVs with pandas (NO R, NO rdata), write
#    the output .rds with `pyreadr`.  The CSVs (simulated_data/sim_data_all.csv,
#    true_params.csv) are produced ONCE from the .rds and committed (see the
#    converter snippet in the deliverables).  Only I/O dependency is `pyreadr`.
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
                          [["day", "dentadult", "lumadult"]]
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


def write_result_rds(result: dict, path: Path) -> None:
    """Write the fit result as an .rds list-equivalent (read by collect_lrt.R).

    collect_lrt.R reads res$ll, res$se, res$coef and the metadata
    res[c("Np","Mp","Np_rep","block")].  pyreadr writes a one-row data.frame
    rather than an R `list`, but `res$ll`, `res$Np` etc. resolve identically on a
    data.frame (column extraction by $), and `res$coef` is recovered from the
    coef_* columns.  See the collect_lrt.R adaptation note in the deliverables.
    """
    import pyreadr
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "ll": [float(result["ll"])],
        "se": [float(result["se"])],
        "Np": [int(result["Np"])],
        "Mp": [int(result["Mp"])],
        "Np_rep": [int(result["Np_rep"])],
        "block": [bool(result["block"])],
        "Nmif_round1": [int(result["Nmif"][0])],
        "Nmif_round2": [int(result["Nmif"][1])],
    }
    # coef as coef_<name> columns; for an alt the unit-specific params are stored
    # per unit as coef_<name>[<unit>] (matching panelPomp::coef naming).
    for k, v in result["coef"].items():
        row[f"coef_{k}"] = [float(v)]
    pyreadr.write_rds(str(path), pd.DataFrame(row))


# ----------------------------------------------------------------------------
# 4. rw_sd  [R: rw_sd(...) fit_null.R lines 144-150]
#    ALL fitted params perturbed at `value`; sigSn = sigSi = 0 always pinned.
# ----------------------------------------------------------------------------
def make_rw_sd(value: float) -> dict[str, float]:
    rw = {name: value for name in PARAM_NAMES}
    rw["sigSn"] = 0.0          # [R: rw_sd(sigSn=0, ...)]
    rw["sigSi"] = 0.0          # [R: rw_sd(sigSi=0, ...)]
    return rw


# ----------------------------------------------------------------------------
# 5. Panel construction + one MPIF round (batched over starts for GPU memory)
# ----------------------------------------------------------------------------
def logmeanexp_se(values: np.ndarray, axis: int):
    """logmeanexp and its Monte-Carlo SE along `axis`  [matches pomp::logmeanexp(se=TRUE)]."""
    m = np.nanmax(values, axis=axis, keepdims=True)
    w = np.exp(values - m)                       # weights in [0,1]
    n = np.sum(~np.isnan(values), axis=axis)
    mean_w = np.nanmean(w, axis=axis)
    est = np.squeeze(m, axis=axis) + np.log(mean_w)
    se = np.nanstd(w, axis=axis, ddof=1) / mean_w / np.sqrt(n)   # delta-method SE
    return est, se


def build_pomp_dict(sim_data: dict[str, pd.DataFrame], theta: dict[str, float]):
    """9 per-unit pp.Pomp objects  [R: pomplist loop  fit_null.R lines 96-118]."""
    pomp_dict = {}
    for unit_name in UNIT_NAMES:
        ys_u = (sim_data[unit_name]
                .set_index("day")[["dentadult", "lumadult"]]
                .astype(float))
        pomp_dict[unit_name] = pp.Pomp(
            ys=ys_u, theta=pp.PompParameters(theta), statenames=SRJF_STATENAMES, t0=1.0,
            rinit=srjf_rinit, rproc=srjf_rproc,
            dmeas=srjf_dmeas, rmeas=srjf_rmeas,
            par_trans=srjf_par_trans, dt=0.25,            # [R: euler(delta.t = 1/4)]
            accumvars=("error_count",),                   # [R: accumvars = c("error_count")]
        )
    return pomp_dict


def make_panel_parameters(param_rows: pd.DataFrame, specific_names: list[str]):
    """
    Build PanelParameters, one spec per start.

      * NULL (specific_names == []):  ALL 15 params shared, EMPTY unit_specific.
        [R: panelPomp(pomplist, shared = true_params)]  -> block=True is a no-op,
        so this is the all-shared PIF.

      * ALT (specific_names != []):  the alt block becomes UNIT-SPECIFIC -- one
        column per unit, all units START at the shared value (warm start = the
        round value), the rest stay shared.
        [R: create_parameters() -> shared = params minus block;
            specific = matrix(block_param, nrow=len(block), ncol=9, byrow=TRUE);
            panelPomp(pomplist, shared=..., specific=...)]
        block=True MPIF then resamples the unit-specific block per unit.
    """
    shared_names = [n for n in PARAM_NAMES if n not in specific_names]
    theta_specs = []
    for _, row in param_rows.iterrows():
        shared_df = pd.DataFrame(
            {"shared": [float(row[name]) for name in shared_names]},
            index=shared_names,
        )
        if specific_names:
            # one row per alt param, one column per unit; all units = the start value
            # [R: matrix(rep(param, 9), nrow=len(block), ncol=9, byrow=TRUE)]
            unit_df = pd.DataFrame(
                {u: [float(row[name]) for name in specific_names] for u in UNIT_NAMES},
                index=specific_names,
            )
        else:
            unit_df = pd.DataFrame(index=[], columns=UNIT_NAMES)
        theta_specs.append({"shared": shared_df, "unit_specific": unit_df})
    return pp.PanelParameters(theta=theta_specs)


def panel_theta_to_coef(theta_dict: dict, specific_names: list[str]) -> dict[str, float]:
    """
    Recover a panelPomp::coef-style flat named vector from one panel-theta spec:
      shared params -> "<name>"
      unit-specific -> "<name>[<unit>]"   (matches panelPomp coef naming)
    """
    coef = {}
    shared_df = theta_dict["shared"]
    for name in shared_df.index:
        coef[name] = float(shared_df.loc[name, "shared"])
    unit_df = theta_dict["unit_specific"]
    if specific_names and len(unit_df.index) > 0:
        for name in unit_df.index:
            for u in UNIT_NAMES:
                coef[f"{name}[{u}]"] = float(unit_df.loc[name, u])
    return coef


def panel_theta_to_shared_row(theta_dict: dict) -> dict[str, float]:
    """Shared params only, for the round-1 -> round-2 warm-start carry.  The
    unit-specific block is NOT collapsed here; it is carried verbatim in the raw
    `theta` spec (see select_round_two) so round 2 resumes the FULL per-unit
    estimates  [R: fit_alt.R:206-211 replicated_specific_list]."""
    sd = theta_dict["shared"]
    return {str(n): float(sd.loc[n, "shared"]) for n in sd.index}


def _run_mif_batch(pomp_dict, starts_batch, specific_names, rw_sd_dict, nmif,
                   mif_particles, pfilter_particles, pfilter_reps,
                   mif_key, pf_key, vmap_chunk_size):
    """MPIF + replicated pfilter on ONE batch of starts -> list of (coef, shared, theta, ll, se)."""
    panel = pp.PanelPomp(
        Pomp_dict=pomp_dict,
        theta=make_panel_parameters(starts_batch, specific_names),
    )

    # [R: mif2(panelfood, block=TRUE, Nmif=nmif, rw.sd=..., cooling.fraction.50=0.7, Np=Mp)]
    # block=True is pypomp's MPIF.  For the NULL (empty unit_specific) it equals PIF
    # (block is a no-op); for an ALT it block-resamples the unit-specific params per
    # unit -- exactly matching the R alt fit (block=TRUE).
    # pypomp >=0.4.6: cooling is carried by RWSigma via .geometric_cooling(a=...),
    # NOT a mif(a=...) kwarg.  a=0.7 == R cooling.fraction.50=0.7.
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

    # logLiks: shape (Nstarts, U, reps) of per-unit, per-rep log-likelihoods.
    # panel_logmeanexp(MARGIN=1) = sum_units logmeanexp_over_reps  (NOT sum-then-lme).
    loglik_array = np.asarray(panel.results_history[-1].logLiks.values)
    assert loglik_array.ndim == 3, f"expected (Nstarts,U,reps), got {loglik_array.shape}"
    per_unit_lme, per_unit_se = logmeanexp_se(loglik_array, axis=2)   # (Nstarts, U) each
    panel_loglik = np.nansum(per_unit_lme, axis=1)                    # (Nstarts,)
    panel_se = np.sqrt(np.nansum(per_unit_se ** 2, axis=1))           # combine unit SEs

    out = []
    for idx, theta_dict in enumerate(panel.theta._to_list()):
        out.append({
            "coef": panel_theta_to_coef(theta_dict, specific_names),   # full coef (output schema)
            "shared": panel_theta_to_shared_row(theta_dict),           # for warm-start carry
            "theta": theta_dict,                                       # raw spec (carries unit-specific)
            "loglik": float(panel_loglik[idx]),
            "se": float(panel_se[idx]),
        })
    return out


def run_mif_round(pomp_dict, starts: pd.DataFrame, specific_names, rw_sd_dict,
                  nmif, mif_particles, pfilter_particles, pfilter_reps,
                  key_seed, vmap_chunk_size):
    """One round over ALL starts, fed to pypomp in batches of BATCH_SIZE."""
    starts = starts.reset_index(drop=True)
    n = len(starts)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    out = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts.iloc[lo:lo + BATCH_SIZE]
        print(f"    batch {bi + 1}/{n_batches}  (starts {lo}..{lo + len(batch) - 1})", flush=True)
        out.extend(_run_mif_batch(
            pomp_dict, batch, specific_names, rw_sd_dict, nmif, mif_particles,
            pfilter_particles, pfilter_reps,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size,
        ))
    return out


def select_round_two(round_one: list, specific_names):
    """
    Top ceil(N/4) by loglik, each replicated x4, carrying BOTH the shared warm
    start AND (for the alt) the per-unit unit-specific block.
    [R: fit_null.R lines 160-167 / fit_alt.R lines 196-211  select top 25% +
    rep(each=4); the per-unit @specific is carried verbatim as
    replicated_specific_list -- NOT collapsed to one value].
    Returns (starts_df, carry list aligned to starts_df rows).
    """
    order = sorted(range(len(round_one)), key=lambda i: round_one[i]["loglik"], reverse=True)
    n_select = math.ceil(len(round_one) / 4.0)
    sel = order[:n_select]
    starts_rows, carry = [], []
    for i in sel:
        shared = round_one[i]["shared"]
        td = round_one[i]["theta"]
        for _ in range(4):                       # replicate x4  [R: replicate(4, ..., simplify=FALSE)]
            # round-2 starts carry the shared params; the unit-specific block is
            # carried separately in `carry` (used to rebuild PanelParameters).
            starts_rows.append(shared)
            carry.append(td)
    starts_df = pd.DataFrame(starts_rows)
    # ensure every PARAM_NAME column exists (shared lacks the alt block names for the
    # alt; fill those from the carried unit-specific u1 value so make_panel_parameters
    # has a full row even though the specific_names columns are then overwritten by carry).
    for n in PARAM_NAMES:
        if n not in starts_df.columns:
            starts_df[n] = [float(carry[r]["unit_specific"].loc[n, UNIT_NAMES[0]])
                            for r in range(len(carry))]
    return starts_df[PARAM_NAMES], carry


def make_panel_parameters_with_carry(starts_df, carry, specific_names):
    """Like make_panel_parameters, but the alt's unit-specific block is taken from
    the carried round-1 per-unit matrices (NOT collapsed to one value), so round 2
    resumes the per-unit estimates  [R: fit_alt.R:224 specific.start = replicated_specific_list]."""
    shared_names = [n for n in PARAM_NAMES if n not in specific_names]
    specs = []
    for r, (_, row) in enumerate(starts_df.iterrows()):
        shared_df = pd.DataFrame({"shared": [float(row[n]) for n in shared_names]},
                                 index=shared_names)
        if specific_names:
            us_prev = carry[r]["unit_specific"]
            unit_df = pd.DataFrame(
                {u: [float(us_prev.loc[p, u]) for p in specific_names] for u in UNIT_NAMES},
                index=specific_names,
            )
        else:
            unit_df = pd.DataFrame(index=[], columns=UNIT_NAMES)
        specs.append({"shared": shared_df, "unit_specific": unit_df})
    return pp.PanelParameters(theta=specs)


def _run_mif_batch_carry(pomp_dict, starts_batch, carry_batch, specific_names, rw_sd_dict,
                         nmif, mif_particles, pfilter_particles, pfilter_reps,
                         mif_key, pf_key, vmap_chunk_size):
    """Round-2 MPIF + replicated pfilter on ONE batch, seeding the alt block per-unit
    from the carried round-1 matrices  [R: fit_alt.R:226-242 round-2 mif2]."""
    panel = pp.PanelPomp(
        Pomp_dict=pomp_dict,
        theta=make_panel_parameters_with_carry(starts_batch, carry_batch, specific_names),
    )
    panel.mif(
        J=mif_particles, M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
        key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size,
    )
    panel.pfilter(
        J=pfilter_particles, reps=pfilter_reps, key=pf_key,
        chunk_size=vmap_chunk_size if vmap_chunk_size else 1,
    )
    loglik_array = np.asarray(panel.results_history[-1].logLiks.values)
    assert loglik_array.ndim == 3, f"expected (Nstarts,U,reps), got {loglik_array.shape}"
    per_unit_lme, per_unit_se = logmeanexp_se(loglik_array, axis=2)
    panel_loglik = np.nansum(per_unit_lme, axis=1)
    panel_se = np.sqrt(np.nansum(per_unit_se ** 2, axis=1))
    out = []
    for idx, theta_dict in enumerate(panel.theta._to_list()):
        out.append({
            "coef": panel_theta_to_coef(theta_dict, specific_names),
            "loglik": float(panel_loglik[idx]),
            "se": float(panel_se[idx]),
        })
    return out


def run_mif_round_carry(pomp_dict, starts_df, carry, specific_names, rw_sd_dict,
                        nmif, mif_particles, pfilter_particles, pfilter_reps,
                        key_seed, vmap_chunk_size):
    """Round-2 driver: like run_mif_round but propagates the per-unit alt-block carry."""
    starts_df = starts_df.reset_index(drop=True)
    n = len(starts_df)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    out = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts_df.iloc[lo:lo + BATCH_SIZE]
        carry_batch = carry[lo:lo + len(batch)]
        print(f"    batch {bi + 1}/{n_batches}  (starts {lo}..{lo + len(batch) - 1})", flush=True)
        out.extend(_run_mif_batch_carry(
            pomp_dict, batch, carry_batch, specific_names, rw_sd_dict, nmif,
            mif_particles, pfilter_particles, pfilter_reps,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size,
        ))
    return out


# ----------------------------------------------------------------------------
# 6. Main: two-round global search  [R: fit_null.R / fit_alt.R]
# ----------------------------------------------------------------------------
def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="GPU/pypomp SRJF2 bootstrap-LRT fit for one (dataset, target).")
    p.add_argument("dataset_index", type=int, help="bootstrap dataset index 1..100")
    p.add_argument("target", type=str,
                   help='"null" (all-shared) or an alternative: '
                        + ", ".join(ALT_PARAMS_MAP.keys()))
    p.add_argument("--run-level", type=int,
                   default=int(os.environ.get("PYPOMP_RUN_LEVEL", "3")),
                   choices=[1, 2, 3], help="particle/iteration budget (default 3 == R)")
    p.add_argument("--n-starts", type=int, default=N_STARTS_DEFAULT,
                   help="number of search starts (R uses 3*100=300)")
    return p.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    b = args.dataset_index
    target = args.target
    rl = args.run_level
    if not (1 <= b <= 100):
        raise ValueError("dataset index must be 1..100")

    is_null = (target == "null")
    if not is_null and target not in ALT_PARAMS_MAP:
        raise ValueError(f"unknown target {target!r}; valid: null, "
                         + ", ".join(ALT_PARAMS_MAP.keys()))
    specific_names = [] if is_null else ALT_PARAMS_MAP[target]

    rl_i = rl - 1
    Mp = ALGORITHMIC_PARAMS["Mp"][rl_i]
    Np = ALGORITHMIC_PARAMS["Np"][rl_i]
    Np_rep = ALGORITHMIC_PARAMS["Np_rep"][rl_i]
    vmap_chunk = os.environ.get("PYPOMP_VMAP_CHUNK")
    vmap_chunk = int(vmap_chunk) if vmap_chunk else None

    print(f"=== pypomp SRJF2 bootstrap LRT: dataset {b}, target {target} "
          f"(run_level {rl}, {args.n_starts} starts) ===")
    print(f"JAX backend = {jax.default_backend()}   x64 = {_USE_X64}   "
          f"devices = {jax.devices()}   vmap_chunk = {vmap_chunk}")
    print(f"Mp(mif J) = {Mp}  Np(pfilter J) = {Np}  Np_rep = {Np_rep}  "
          f"batch = {BATCH_SIZE}  unit_specific = {specific_names or 'none (null)'}")

    base = Path(__file__).resolve().parent
    sim_dir = base / "simulated_data"
    # Output path  [R: results_null/lrt_null_<b>.rds | results_alt/lrt_<alt>_<b>.rds].
    # PYPOMP_OVERWRITE=0 (default) writes a *_pypomp.rds SIDECAR (does NOT touch the
    # canonical R result); PYPOMP_OVERWRITE=1 overwrites the canonical lrt_*.rds.
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

    # Round-1 starts: all = the true params (warm start)  [R: shared.start = true_params].
    starts = pd.DataFrame([theta0] * args.n_starts)[PARAM_NAMES]

    # seed offset distinguishes null vs each alt  [R: set.seed(801*1000 + b [+ 100000 for alt])]
    target_idx = 0 if is_null else (list(ALT_PARAMS_MAP.keys()).index(target) + 1)
    seed_base = 801_000 + b + 1000 * target_idx

    print("Round 1 starting...", flush=True)
    round_one = run_mif_round(
        pomp_dict, starts, specific_names, make_rw_sd(RW_ROUND1),
        nmif=NMIF_ROUND1, mif_particles=Mp, pfilter_particles=Np,
        pfilter_reps=Np_rep, key_seed=seed_base + 10_000,
        vmap_chunk_size=vmap_chunk,
    )
    print(f"Round 1 done.  best ll = {max(r['loglik'] for r in round_one):.3f}", flush=True)

    # Select top 25%, replicate x4, carrying the FULL per-unit unit-specific block
    # [R: fit_alt.R:196-211 replicated_specific_list].
    round_two_starts, carry = select_round_two(round_one, specific_names)

    print("Round 2 starting...", flush=True)
    round_two = run_mif_round_carry(
        pomp_dict, round_two_starts, carry, specific_names, make_rw_sd(RW_ROUND2),
        nmif=NMIF_ROUND2, mif_particles=Mp, pfilter_particles=Np,
        pfilter_reps=Np_rep, key_seed=seed_base + 20_000,
        vmap_chunk_size=vmap_chunk,
    )

    # [R: lls <- ...; best <- which.max(lls[1,]); result <- list(ll, se, coef, ...)]
    best = max(round_two, key=lambda r: r["loglik"])
    result = {
        "ll": best["loglik"],
        "se": best["se"],
        "coef": best["coef"],
        "Np": Np,
        "Mp": Mp,
        "Np_rep": Np_rep,
        "block": True,
        "Nmif": (NMIF_ROUND1, NMIF_ROUND2),    # round1=150, round2=250
    }
    write_result_rds(result, out_path)

    print("\n================= RESULT =================")
    print(f"target = {target}   dataset b = {b}")
    print(f"ll = {result['ll']:.4f}   se = {result['se']:.4f}")
    print(f"saved -> {out_path}")
    print("==========================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
