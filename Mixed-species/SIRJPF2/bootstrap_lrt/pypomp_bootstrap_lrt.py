#!/usr/bin/env python3
"""
pypomp_bootstrap_lrt.py
=======================
GPU / pypomp port of the SIRJPF2 parametric-bootstrap-LRT computation
(fit_null.R + fit_alt.R).  ONE script that fits, for a given bootstrap
replicate ``b`` and a given ``target`` (the all-shared NULL, or one of the
seven ALTERNATIVES whose parameter block is made unit-specific), the same
two-round marginalized panel iterated filtering (MPIF) fit on the SAME
simulated panel that the R scripts use, then pfilters and writes the result
in the EXACT schema that collect_lrt.R consumes.

This is the REFERENCE family for the bootstrap-LRT port and doubles as the
VALIDATION GATE: the SIRJPF2 JAX model is already validated in
``../model/pypomp_all_shared.py`` (it reproduced the R all-shared MLE ll,
-880.44 vs R -880.56, on GPU).  This script REUSES that model verbatim and
adds only the bootstrap-LRT driver (warm starts, the alt unit-specific
block, and the collect_lrt.R output schema).  See "VALIDATION" below for
how to check the per-b ll against the existing R results.

Usage
-----
    python pypomp_bootstrap_lrt.py <b 1..100> <target>

    target = "null"                 -> all-shared fit  -> results_null/lrt_null_<b>.rds
           | one of the 7 alts      -> block fit       -> results_alt/lrt_<alt>_<b>.rds

    alternatives (block becomes unit-specific):
        xi                  -> xi
        theta_Si_theta_Sn   -> theta_Si, theta_Sn
        theta_P             -> theta_P
        theta_Ii_theta_In   -> theta_Ii, theta_In
        ri_rn               -> ri, rn
        probi_probn         -> probi, probn
        f_Si_f_Sn           -> f_Si, f_Sn

    examples
        python pypomp_bootstrap_lrt.py 7 null
        python pypomp_bootstrap_lrt.py 7 xi
        python pypomp_bootstrap_lrt.py 42 ri_rn

Environment toggles (defaults reproduce the R run):
    PYPOMP_USE_GPU   = 1 (default) | 0/cpu     -> JAX device
    PYPOMP_X64       = 1 (default) | 0          -> float64 (R uses double; keep 1)
    PYPOMP_RUN_LEVEL = 3 (default) | 1 | 2      -> particle/iteration budget
                                                   (1 = tiny CPU smoke ONLY)
    PYPOMP_NSTARTS   = 300 (default)            -> number of mif starts
                                                   (R uses 3*getDoParWorkers()=3*100=300)
    PYPOMP_BATCH     = 100 (default) | <int>    -> outer batch over starts (GPU memory)
    PYPOMP_VMAP_CHUNK= (unset) | <int>          -> pypomp per-UNIT vmap chunk

REQUIRES the development pypomp (>=0.4.6, with PanelPomp/PanelParameters/
RWSigma/ParTrans) and a CUDA JAX on the GPU host.  Inputs are the committed
CSVs (../coverage_study/simulated_data/sim_data_all.csv + true_params.csv,
identical to the .rds the R scripts read), read with pandas, so NO R and NO
rdata reader is needed on the GPU node.  The only output-I/O dependency is
``pyreadr`` (writes the .rds in collect_lrt.R's schema):  pip install pyreadr.

================================  FIDELITY NOTE  ===============================
JAX's PRNG is NOT R's, so this is NOT bit-identical: the particle draws and the
mif perturbations differ from a run of fit_null.R / fit_alt.R.  What this
reproduces EXACTLY is the estimator + workflow:
  * identical SIRJPF2 process / measurement model (term-for-term with the R
    Csnippets; the model functions are byte-identical to pypomp_all_shared.py,
    itself GPU-validated against the R MLE),
  * the same 2-round MPIF: Nmif 150 then 250, Np=Mp=1500, Np_rep=10,
    cooling.fraction.50=0.7, top-25%-and-replicate-x4 selection, block=True,
  * the alt's parameter block made UNIT-SPECIFIC (non-empty PanelParameters
    unit_specific) exactly as create_parameters() does in fit_alt.R,
  * the panel loglik = sum_units logmeanexp_over_reps(unitLogLik) with its SE,
    matching panelPomp::panel_logmeanexp(..., MARGIN=1, se=TRUE).
Compare the two STATISTICALLY (per-b ll within the pfilter MC SE), NOT
byte-for-byte.  R<->pypomp mapping annotated inline as  [R: <file>:<line>].

================================  VALIDATION  =================================
The gate: for a few b, the pypomp ll must match the R ll within the pfilter
Monte-Carlo SE.  Procedure (CPU or GPU, run_level 3):

    # pick e.g. b in {1, 7, 42}
    for b in 1 7 42; do
        python pypomp_bootstrap_lrt.py $b null
        python pypomp_bootstrap_lrt.py $b xi
    done

Then in R compare against the EXISTING R fits:

    for (b in c(1,7,42)) {
      rn <- readRDS(sprintf("results_null/lrt_null_%d.rds", b))      # R null
      pn <- readRDS(sprintf("results_null/lrt_null_%d_pypomp.rds",b))# pypomp null
      cat(sprintf("b=%d NULL  R ll=%.2f (se %.2f)  pypomp ll=%.2f (se %.2f)  d=%.2f\n",
          b, rn$ll, rn$se, pn$ll, pn$se, pn$ll-rn$ll))
      ra <- readRDS(sprintf("results_alt/lrt_xi_%d.rds", b))
      pa <- readRDS(sprintf("results_alt/lrt_xi_%d_pypomp.rds", b))
      cat(sprintf("b=%d ALT   R ll=%.2f (se %.2f)  pypomp ll=%.2f (se %.2f)  d=%.2f\n",
          b, ra$ll, ra$se, pa$ll, pa$se, pa$ll-ra$ll))
    }

PASS criterion: |pypomp ll - R ll| within ~2-3x the combined pfilter MC SE for
both arms (the ll is a stochastic MLE, so small gaps are expected; both are
estimating the same maximized panel loglik).  The bootstrap Lambda = 2*(ll_alt
- ll_null), so the GATE is really that BOTH arms land in the same place as R.
By default pypomp writes a *_pypomp.rds sidecar (set PYPOMP_OVERWRITE=1 to
write the canonical lrt_*.rds name and feed collect_lrt.R directly).
===============================================================================
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

_USE_GPU = os.environ.get("PYPOMP_USE_GPU", "1").lower() not in {"0", "false", "no", "cpu"}
if _USE_GPU:
    os.environ.pop("JAX_PLATFORMS", None)
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
N_UNITS = 8
UNIT_NAMES = [f"u{i}" for i in range(1, N_UNITS + 1)]   # [R: fit_null.R:116 1:8 -> u1..u8]

# Canonical 26-parameter order  [R: fit_null.R:106-109 param_names]
PARAM_NAMES = [
    "xi", "sigSn", "sigIn", "sigSi", "sigIi", "sigF", "sigP",
    "theta_Sn", "theta_In", "theta_Si", "theta_P", "theta_Ii",
    "f_Sn", "f_Si", "rn", "ri", "probn", "probi",
    "k_Ii", "k_In", "k_Sn", "k_Si", "sigJi", "sigJn", "theta_Jn", "theta_Ji",
]

# State vector  [R: fit_null.R:110 state_names];  order irrelevant (dict-keyed).
SIRJPF_STATENAMES = [
    "Sn", "In", "Jn", "Si", "Ii", "Ji", "F", "P",
    "T_Sn", "T_In", "T_Si", "T_Ii", "error_count",
]

# Log-transformed params  [R: fit_null.R:100-104 parameter_trans(log=...)].
# sigSn, sigSi are fixed at 0 (rw.sd = 0) -- NOT log-transformed (log(0) = -inf);
# in R they ARE in the log= list but rw.sd=0 keeps them pinned at the start 0.
SIRJPF_LOG_PARAMS = (
    "rn", "ri", "f_Sn", "f_Si", "probn", "probi", "xi",
    "theta_Sn", "theta_Si", "theta_In", "theta_Ii", "theta_Jn", "theta_Ji", "theta_P",
    "sigIn", "sigIi", "sigJn", "sigJi", "sigF", "sigP", "k_Sn", "k_Si", "k_In", "k_Ii",
)
FIXED_ZERO_PARAMS = ("sigSn", "sigSi")   # [R: rw_sd(sigSn=0, sigSi=0); pinned at 0]

# ---- Algorithmic params  [R: fit_null.R:143-145 / fit_alt.R:180-182] ----
# Np=Mp=1500, Np_rep=10 (raised from 1000), block=TRUE for alts, cooling a=0.7.
# Nmif round1=150, round2=250 (hard-coded in the mif2 calls).
NP_FINAL = 1500          # [R: Np = 1500]  pfilter particles
NP_REP = 10              # [R: Np_rep = 10]  pfilter replicates
MP = 1500                # [R: Mp = 1500]  mif2 particles
NMIF_ROUNDS = (150, 250)  # [R: Nmif = 150 (round1), 250 (round2)]
RW_ROUNDS = (0.05, 0.04)  # [R: dent_rw.sd = 0.05 then 0.04]
COOLING_A = 0.7          # [R: cooling.fraction.50 = 0.7]
N_STARTS_DEFAULT = int(os.environ.get("PYPOMP_NSTARTS", "300"))  # [R: 3*getDoParWorkers()=3*100]
BATCH_SIZE = int(os.environ.get("PYPOMP_BATCH", "100"))

# Run-level budgets (run_level 1 = tiny CPU smoke ONLY; 3 = the real R settings).
ALGO = {
    "Np":     {1: 50, 2: 500, 3: NP_FINAL},
    "Np_rep": {1: 2,  2: 10,  3: NP_REP},
    "Mp":     {1: 50, 2: 500, 3: MP},
    "Nmif":   {1: (3, 4), 2: (60, 100), 3: NMIF_ROUNDS},
    "nstarts": {1: 4, 2: 60, 3: N_STARTS_DEFAULT},
}

# ---- The 7 alternatives: name -> the block of params made UNIT-SPECIFIC ----
# [R: fit_alt.R:18-26 alt_params_map]
ALT_PARAMS_MAP = {
    "xi":                ["xi"],
    "theta_Si_theta_Sn": ["theta_Si", "theta_Sn"],
    "theta_P":           ["theta_P"],
    "theta_Ii_theta_In": ["theta_Ii", "theta_In"],
    "ri_rn":             ["ri", "rn"],
    "probi_probn":       ["probi", "probn"],
    "f_Si_f_Sn":         ["f_Si", "f_Sn"],
}
# seed offset per alt (stable ordering for the PRNG keys)
ALT_ORDER = list(ALT_PARAMS_MAP.keys())


# ============================================================================
# 2. SIRJPF2 model in JAX  -- BYTE-IDENTICAL to ../model/pypomp_all_shared.py,
#    which is term-for-term with the fit_null.R / fit_alt.R Csnippets and was
#    GPU-validated against the R all-shared MLE (-880.44 vs R -880.56).
#    R<->JAX term mapping below cites fit_null.R lines (fit_alt.R is identical).
# ============================================================================
def sirjpf_rproc(X_, theta_, key, covars, t, dt):
    """One Euler-Maruyama step  [R: fit_null.R:21-68 euler(delta.t = 1/4)]."""
    Sn, Jn, In = X_["Sn"], X_["Jn"], X_["In"]
    Si, Ji, Ii = X_["Si"], X_["Ji"], X_["Ii"]
    F, P = X_["F"], X_["P"]
    error_count = X_["error_count"]
    sigSn, sigIn = theta_["sigSn"], theta_["sigIn"]
    sigSi, sigIi = theta_["sigSi"], theta_["sigIi"]
    sigJn, sigJi = theta_["sigJn"], theta_["sigJi"]
    sigF, sigP = theta_["sigF"], theta_["sigP"]
    theta_Sn, theta_In = theta_["theta_Sn"], theta_["theta_In"]
    theta_Si, theta_Ii = theta_["theta_Si"], theta_["theta_Ii"]
    theta_Jn, theta_Ji = theta_["theta_Jn"], theta_["theta_Ji"]
    theta_P = theta_["theta_P"]
    f_Sn, f_Si = theta_["f_Sn"], theta_["f_Si"]
    rn, ri = theta_["rn"], theta_["ri"]
    probn, probi = theta_["probn"], theta_["probi"]
    xi = theta_["xi"]

    delta, mu_food, lambda_J = 0.013, 0.37, 0.1   # [R: delta=0.013; +0.37*dt; 0.1*J maturation]
    keys = jax.random.split(key, 8)
    sqdt = jnp.sqrt(dt)
    # [R:26-33  noiX = rnorm(0, sigX * sqrt(dt))]
    noiSn = sigSn * sqdt * jax.random.normal(keys[0])
    noiIn = sigIn * sqdt * jax.random.normal(keys[1])
    noiSi = sigSi * sqdt * jax.random.normal(keys[2])
    noiIi = sigIi * sqdt * jax.random.normal(keys[3])
    noiJn = sigJn * sqdt * jax.random.normal(keys[4])
    noiJi = sigJi * sqdt * jax.random.normal(keys[5])
    noiF = sigF * sqdt * jax.random.normal(keys[6])
    noiP = sigP * sqdt * jax.random.normal(keys[7])

    # [R:35-42  the eight *_term drift+diffusion increments]
    Sn_term = (lambda_J * Jn * dt - theta_Sn * Sn * dt
               - probn * f_Sn * Sn * P * dt - delta * Sn * dt + Sn * noiSn)
    Jn_term = (rn * f_Sn * F * Sn * dt - lambda_J * Jn * dt
               - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn)
    In_term = (probn * f_Sn * Sn * P * dt - theta_In * In * dt
               - delta * In * dt + In * noiIn)
    Si_term = (lambda_J * Ji * dt - theta_Si * Si * dt
               - probi * f_Si * Si * P * dt - delta * Si * dt + Si * noiSi)
    Ji_term = (ri * f_Si * F * Si * dt - lambda_J * Ji * dt
               - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi)
    Ii_term = (probi * f_Si * Si * P * dt - theta_Ii * Ii * dt
               - delta * Ii * dt + Ii * noiIi)
    F_term = (-f_Sn * F * (Sn + xi * In + Jn) * dt
              - f_Si * F * (Si + xi * Ii + Ji) * dt
              - delta * F * dt + mu_food * dt + F * noiF)
    P_term = (30.0 * theta_In * In * dt + 30.0 * theta_Ii * Ii * dt
              - f_Sn * (Sn + xi * In) * P * dt - f_Si * (Si + xi * Ii) * P * dt
              - theta_P * P * dt - delta * P * dt + P * noiP)

    Sn_new = Sn + Sn_term; In_new = In + In_term; Jn_new = Jn + Jn_term
    Si_new = Si + Si_term; Ii_new = Ii + Ii_term; Ji_new = Ji + Ji_term
    F_new = F + F_term; P_new = P + P_term
    # [R:53  if (t-4 ~ 0) P += 25]  -- parasite inoculation pulse at day 4
    P_new = P_new + jnp.where(jnp.abs(t - 4.0) < 0.001, 25.0, 0.0)

    # [R:55-62  soft constraints -> error_count increments, then state clamped to 0]
    def viol(x, hi):
        return (x < 0.0) | (x > hi)

    eps = 0.0
    eps += jnp.where(viol(Sn_new, 1e5), 1.0, 0.0)        # [R:55 += 1]
    eps += jnp.where(viol(Si_new, 1e5), 1.0e6, 0.0)      # [R:56 += 1000000]
    eps += jnp.where(viol(F_new, 1e20), 1.0e3, 0.0)      # [R:57 += 1000]
    eps += jnp.where(viol(In_new, 1e5), 1.0e-3, 0.0)     # [R:58 += 0.001]
    eps += jnp.where(viol(Ii_new, 1e5), 1.0e-9, 0.0)     # [R:59 += 1e-9]
    eps += jnp.where(viol(Jn_new, 1e5), 1.0e-3, 0.0)     # [R:60 += 0.001]
    eps += jnp.where(viol(Ji_new, 1e5), 1.0e-9, 0.0)     # [R:61 += 1e-9]
    eps += jnp.where(viol(P_new, 1e20) & (t > 3.9), 1.0e-6, 0.0)  # [R:62 += 1e-6 if t>3.9]

    Sn_new = jnp.clip(Sn_new, 0.0, 1e5); In_new = jnp.clip(In_new, 0.0, 1e5)
    Jn_new = jnp.clip(Jn_new, 0.0, 1e5); Si_new = jnp.clip(Si_new, 0.0, 1e5)
    Ii_new = jnp.clip(Ii_new, 0.0, 1e5); Ji_new = jnp.clip(Ji_new, 0.0, 1e5)
    F_new = jnp.clip(F_new, 0.0, 1e20); P_new = jnp.clip(P_new, 0.0, 1e20)

    return {
        "Sn": Sn_new, "In": In_new, "Jn": Jn_new, "Si": Si_new, "Ii": Ii_new,
        "Ji": Ji_new, "F": F_new, "P": P_new,
        # [R:64-67  T_Sn = fabs(Sn), ...]  measurement-ready non-negative densities
        "T_Sn": jnp.abs(Sn_new), "T_In": jnp.abs(In_new),
        "T_Si": jnp.abs(Si_new), "T_Ii": jnp.abs(Ii_new),
        # [R: error_count is an accumvar -> reset to 0 each obs by pypomp accumvars]
        "error_count": error_count + eps,
    }


def sirjpf_rinit(theta_, key, covars, t0):
    """Deterministic initial state  [R: fit_null.R:70-79 dyn_init]."""
    return {
        "Sn": jnp.array(2.333), "In": jnp.array(0.0), "Jn": jnp.array(0.0),
        "Si": jnp.array(0.667), "Ii": jnp.array(0.0), "Ji": jnp.array(0.0),
        "F": jnp.array(16.667), "P": jnp.array(0.0),
        "T_Sn": jnp.array(0.0), "T_In": jnp.array(0.0),
        "T_Si": jnp.array(0.0), "T_Ii": jnp.array(0.0),
        "error_count": jnp.array(0.0),
    }


def _nb_logpmf(y, mu, size):
    """Negative-binomial log-pmf in (mu, size) form  [R: dnbinom_mu(y, size, mu, log)]."""
    mu = jnp.maximum(mu, 1e-10); size = jnp.maximum(size, 1e-10)
    return (gammaln(y + size) - gammaln(size) - gammaln(y + 1.0)
            + size * jnp.log(size / (size + mu)) + y * jnp.log(mu / (size + mu)))


def sirjpf_dmeas(Y_, X_, theta_, covars, t):
    """4-D NB measurement log-likelihood  [R: fit_null.R:81-91 dmeas]."""
    ll = (_nb_logpmf(Y_["dentadult"], X_["T_Sn"], theta_["k_Sn"])
          + _nb_logpmf(Y_["dentinf"], X_["T_In"], theta_["k_In"])
          + _nb_logpmf(Y_["lumadult"], X_["T_Si"], theta_["k_Si"])
          + _nb_logpmf(Y_["luminf"], X_["T_Ii"], theta_["k_Ii"]))
    # [R:82-83  if (error_count > 0.0) lik = -150]
    return jnp.where(X_["error_count"] > 0.0, -150.0, ll)


def _nb_sample(key, mu, size):
    mu = jnp.maximum(mu, 1e-10); size = jnp.maximum(size, 1e-10)
    k1, k2 = jax.random.split(key)
    return jax.random.poisson(k2, jax.random.gamma(k1, size) * (mu / size))


def sirjpf_rmeas(X_, theta_, key, covars, t):
    """Measurement simulator  [R: fit_null.R:93-98 rmeas]  (unused by the LRT fit)."""
    keys = jax.random.split(key, 4)
    return jnp.array([
        _nb_sample(keys[0], X_["T_Sn"], theta_["k_Sn"]),
        _nb_sample(keys[1], X_["T_In"], theta_["k_In"]),
        _nb_sample(keys[2], X_["T_Si"], theta_["k_Si"]),
        _nb_sample(keys[3], X_["T_Ii"], theta_["k_Ii"]),
    ], dtype=float)


def sirjpf_to_est(theta):
    """natural -> estimation scale  [R: parameter_trans(log=...) toEst]."""
    out = {**theta}
    for n in SIRJPF_LOG_PARAMS:
        out[n] = jnp.log(jnp.maximum(theta[n], 1e-30))
    for n in FIXED_ZERO_PARAMS:
        out[n] = theta[n]   # pinned at 0, not log-transformed
    return out


def sirjpf_from_est(theta):
    """estimation -> natural scale  [R: parameter_trans fromEst]."""
    out = {**theta}
    for n in SIRJPF_LOG_PARAMS:
        out[n] = jnp.exp(theta[n])
    for n in FIXED_ZERO_PARAMS:
        out[n] = theta[n]
    return out


sirjpf_par_trans = pp.ParTrans(to_est=sirjpf_to_est, from_est=sirjpf_from_est)


# ============================================================================
# 3. R-FREE I/O: read the SAME inputs the R scripts read, but from the
#    pre-converted CSVs in ../coverage_study/simulated_data (committed; the
#    .rds the GPU node cannot read).  Output via pyreadr in collect_lrt.R's
#    schema.  See make_sim_data_all_csv() below if the CSVs are ever missing.
# ============================================================================
def _sim_dir() -> Path:
    """Directory holding sim_data_all.csv + true_params.csv (shared w/ coverage)."""
    base = Path(__file__).resolve().parent
    d = base.parent / "coverage_study" / "simulated_data"   # [R: ../coverage_study/simulated_data]
    return d


def read_sim_data(b: int) -> dict[str, pd.DataFrame]:
    """Read the b-th simulated panel  [R: fit_null.R:17 readRDS(sim_data_<b>.rds)].
    The CSV (cols: b,unit,day,4 obs) is the EXACT content of those .rds files,
    pre-converted once by the coverage port (see make_sim_data_all_csv)."""
    csv = _sim_dir() / "sim_data_all.csv"
    df = pd.read_csv(csv)
    df = df[df["b"] == b]
    if df.empty:
        raise ValueError(f"no rows for dataset b={b} in {csv}")
    out = {}
    for u in UNIT_NAMES:
        out[u] = (df[df["unit"] == u]
                  [["day", "dentadult", "dentinf", "lumadult", "luminf"]]
                  .astype(float).sort_values("day").reset_index(drop=True))
    return out


def read_true_params() -> dict[str, float]:
    """Read true_params (name,value)  [R: fit_null.R:18 readRDS(true_params.rds)]."""
    tp = pd.read_csv(_sim_dir() / "true_params.csv")
    params = dict(zip(tp["name"].astype(str), tp["value"].astype(float)))
    missing = [n for n in PARAM_NAMES if n not in params]
    if missing:
        raise ValueError(f"true_params.csv missing parameters: {missing}")
    return params


def make_sim_data_all_csv() -> None:
    """ONE-TIME converter (run on a host WITH R) producing the two CSVs this
    script reads from the .rds the R scripts use.  Already done for SIRJPF2
    (the files exist), so this is here only for provenance / regeneration.

    Run with R available, from this directory:

        Rscript --vanilla -e '
          d <- "../coverage_study/simulated_data"
          tp <- readRDS(file.path(d,"true_params.rds"))
          write.csv(data.frame(name=names(tp), value=as.numeric(tp)),
                    file.path(d,"true_params.csv"), row.names=FALSE)
          rows <- list()
          for (b in 1:100) {
            f <- sprintf("%s/sim_data_%03d.rds", d, b)
            if (!file.exists(f)) next
            sl <- readRDS(f)
            for (i in 1:8) {
              u <- sprintf("u%d", i); dat <- as.data.frame(sl[[u]])
              colnames(dat) <- c("day","dentadult","dentinf","lumadult","luminf")
              dat$b <- b; dat$unit <- u
              rows[[length(rows)+1]] <- dat[,c("b","unit","day",
                  "dentadult","dentinf","lumadult","luminf")]
            }
          }
          write.csv(do.call(rbind, rows),
                    file.path(d,"sim_data_all.csv"), row.names=FALSE)'
    """
    raise NotImplementedError("CSVs already exist; see docstring to regenerate with R.")


# ============================================================================
# 4. rw_sd  [R: fit_null.R:163-169 rw_sd(...)]  -- ALL params perturbed except
#    sigSn/sigSi (fixed at 0).  Same set for null and alt; for the alt, the
#    unit-specific block params share this same rw_sd (applied per unit).
# ============================================================================
def make_rw_sd(value: float) -> dict[str, float]:
    rw = {name: value for name in PARAM_NAMES}
    rw["sigSn"] = 0.0   # [R: rw_sd(sigSn=0)]
    rw["sigSi"] = 0.0   # [R: rw_sd(sigSi=0)]
    return rw


# ============================================================================
# 5. Panel construction + the alt's unit-specific block
# ============================================================================
def logmeanexp_se(values: np.ndarray, axis: int):
    """logmeanexp and its delta-method MC SE along `axis`  [R: panel_logmeanexp se=TRUE]."""
    m = np.nanmax(values, axis=axis, keepdims=True)
    w = np.exp(values - m)
    n = np.sum(~np.isnan(values), axis=axis)
    mean_w = np.nanmean(w, axis=axis)
    est = np.squeeze(m, axis=axis) + np.log(mean_w)
    se = np.nanstd(w, axis=axis, ddof=1) / mean_w / np.sqrt(n)
    return est, se


def build_pomp_dict(sim_data, theta):
    """8 per-unit pp.Pomp objects  [R: fit_null.R:115-137 pomplist loop]."""
    pomp_dict = {}
    for u in UNIT_NAMES:
        ys_u = (sim_data[u].set_index("day")[["dentadult", "dentinf", "lumadult", "luminf"]]
                .astype(float))
        pomp_dict[u] = pp.Pomp(
            ys=ys_u, theta=pp.PompParameters(theta), statenames=SIRJPF_STATENAMES, t0=1.0,
            rinit=sirjpf_rinit, rproc=sirjpf_rproc, dmeas=sirjpf_dmeas,
            rmeas=sirjpf_rmeas, par_trans=sirjpf_par_trans, dt=0.25,   # [R: euler delta.t=1/4]
            accumvars=("error_count",),    # [R: accumvars=c("error_count")]
        )
    return pomp_dict


def make_panel_parameters(param_rows: pd.DataFrame, specific_names: list[str]):
    """
    Build a PanelParameters with one spec per mif start.

    NULL  (target == "null", specific_names == []):
        empty unit_specific; all 26 params shared.  block=True is then a no-op
        (PIF == MPIF without unit-specific params), matching fit_null.R.

    ALT  (specific_names == the alt block, e.g. ["xi"] or ["ri","rn"]):
        shared  = the 26 params MINUS the block  [R: create_parameters(): shared =
                  parameters[!(names %in% parameter_names)]],
        unit_specific = the block params, ONE COLUMN PER UNIT (u1..u8), each
                  column initialized to the SAME warm-start value (per-row).
                  [R: fit_alt.R:37-48 create_parameters(): specific_mat is
                  rep(parameters[[p]], 8) -> identical across units at the start,
                  then mif2(block=TRUE) resamples them per unit].
    Each input row is one independent mif start (warm start), vmapped on GPU.
    """
    shared_names = [n for n in PARAM_NAMES if n not in specific_names]
    specs = []
    for _, row in param_rows.iterrows():
        shared_df = pd.DataFrame({"shared": [float(row[n]) for n in shared_names]},
                                 index=shared_names)
        if specific_names:
            # one column per unit; warm start = same value across units (R rep(.,8))
            us = pd.DataFrame(
                {u: [float(row[p]) for p in specific_names] for u in UNIT_NAMES},
                index=specific_names,
            )
        else:
            us = pd.DataFrame(index=[], columns=UNIT_NAMES)
        specs.append({"shared": shared_df, "unit_specific": us})
    return pp.PanelParameters(theta=specs)


def panel_theta_to_coef(theta_dict, specific_names: list[str]) -> dict[str, float]:
    """
    Extract the fitted coefficients in collect_lrt.R's `coef()` schema:
      * shared params by bare name,
      * unit-specific block params as "<param>[u1]".."<param>[u8]".
    [R: panelPomp::coef() naming, e.g. "xi[u1]".."xi[u8]" + shared rest].
    """
    out = {}
    shared_df = theta_dict["shared"]
    for n in shared_df.index:
        out[str(n)] = float(shared_df.loc[n, "shared"])
    us = theta_dict.get("unit_specific")
    if specific_names and us is not None and len(us.index) > 0:
        for p in specific_names:
            for u in UNIT_NAMES:
                out[f"{p}[{u}]"] = float(us.loc[p, u])
    return out


def panel_theta_to_shared_row(theta_dict) -> dict[str, float]:
    """Shared params only, for the round-1 -> round-2 warm-start carry."""
    sd = theta_dict["shared"]
    return {str(n): float(sd.loc[n, "shared"]) for n in sd.index}


# ============================================================================
# 6. One MPIF round (batched over starts) + selection
# ============================================================================
def _run_mif_batch(pomp_dict, starts_batch, specific_names, rw_sd_dict, nmif,
                   Mp, Np, Np_rep, mif_key, pf_key, vmap_chunk_size):
    panel = pp.PanelPomp(
        Pomp_dict=pomp_dict,
        theta=make_panel_parameters(starts_batch, specific_names),
    )
    # [R: mif2(Nmif, rw.sd, cooling.fraction.50=0.7, Np=Mp, block=TRUE)]
    # block=True is pypomp's MPIF.  For the null (empty unit_specific) it equals
    # PIF (block is a no-op); for the alt it block-resamples the unit-specific
    # block -- exactly the fit_alt.R "block=TRUE" branch.
    # pypomp >=0.4.6: cooling rides on RWSigma.geometric_cooling(a=...);
    # a=0.7 == R cooling.fraction.50=0.7.
    panel.mif(J=Mp, M=nmif,
              rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
              key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size)
    # [R: ll <- replicate(Np_rep, unitLogLik(pfilter(m1, Np))); panel_logmeanexp(MARGIN=1, se=TRUE)]
    panel.pfilter(J=Np, reps=Np_rep, key=pf_key,
                  chunk_size=vmap_chunk_size if vmap_chunk_size else 1)
    ll = np.asarray(panel.results_history[-1].logLiks.values)   # (Nstarts, U, reps)
    assert ll.ndim == 3, f"expected (Nstarts,U,reps), got {ll.shape}"
    # logmeanexp over reps PER UNIT, THEN sum over units (does NOT commute with
    # sum-then-logmeanexp).  [R: panel_logmeanexp(x=ll, MARGIN=1, se=TRUE)]
    unit_est, unit_se = logmeanexp_se(ll, axis=2)               # (Nstarts, U) each
    panel_ll = np.nansum(unit_est, axis=1)                      # sum over units
    panel_se = np.sqrt(np.nansum(unit_se ** 2, axis=1))         # combine unit SEs
    rows = []
    theta_list = panel.theta._to_list()
    for idx, td in enumerate(theta_list):
        rows.append({
            "loglik": float(panel_ll[idx]),
            "se": float(panel_se[idx]),
            "coef": panel_theta_to_coef(td, specific_names),     # full coef (output schema)
            "shared": panel_theta_to_shared_row(td),             # for warm-start carry
            "theta": td,                                         # raw spec (carries unit-specific)
        })
    return rows


def run_mif_round(pomp_dict, starts, specific_names, rw_sd_dict, nmif,
                  Mp, Np, Np_rep, key_seed, vmap_chunk_size):
    starts = starts.reset_index(drop=True)
    n = len(starts)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    rows = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts.iloc[lo:lo + BATCH_SIZE]
        print(f"    batch {bi + 1}/{n_batches} (starts {lo}..{lo + len(batch) - 1})", flush=True)
        rows.extend(_run_mif_batch(
            pomp_dict, batch, specific_names, rw_sd_dict, nmif, Mp, Np, Np_rep,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size))
    return rows


def select_round_two(round_one, specific_names):
    """
    Top ceil(N/4) by loglik, each replicated x4, carrying BOTH the shared warm
    start AND (for the alt) the per-unit unit-specific block.
    [R: fit_null.R:179-186 / fit_alt.R:218-233  select top 25% + rep(each=4)].
    Returns (starts_df, specific_carry list aligned to starts_df rows).
    """
    order = sorted(range(len(round_one)), key=lambda i: round_one[i]["loglik"], reverse=True)
    n_select = math.ceil(len(round_one) / 4.0)
    sel = order[:n_select]
    starts_rows, carry = [], []
    for i in sel:
        shared = round_one[i]["shared"]
        td = round_one[i]["theta"]
        for _ in range(4):                       # replicate x4
            # round-2 starts carry the shared params; the unit-specific block is
            # carried separately in `carry` (used to rebuild PanelParameters).
            starts_rows.append(shared)
            carry.append(td)
    starts_df = pd.DataFrame(starts_rows)
    # ensure every PARAM_NAME column exists (shared lacks the alt block names for the
    # alt; fill those from the carried unit-specific u1 value so make_panel_parameters
    # has a full row even though specific_names columns are then overwritten by carry).
    for n in PARAM_NAMES:
        if n not in starts_df.columns:
            starts_df[n] = [float(carry[r]["unit_specific"].loc[n, UNIT_NAMES[0]])
                            for r in range(len(carry))]
    return starts_df[PARAM_NAMES], carry


def make_panel_parameters_with_carry(starts_df, carry, specific_names):
    """Like make_panel_parameters, but the alt's unit-specific block is taken
    from the carried round-1 per-unit matrices (NOT collapsed to one value),
    so round-2 resumes the per-unit estimates  [R: fit_alt.R:246 specific.start]."""
    shared_names = [n for n in PARAM_NAMES if n not in specific_names]
    specs = []
    for r, (_, row) in enumerate(starts_df.iterrows()):
        shared_df = pd.DataFrame({"shared": [float(row[n]) for n in shared_names]},
                                 index=shared_names)
        if specific_names:
            us_prev = carry[r]["unit_specific"]
            us = pd.DataFrame(
                {u: [float(us_prev.loc[p, u]) for p in specific_names] for u in UNIT_NAMES},
                index=specific_names,
            )
        else:
            us = pd.DataFrame(index=[], columns=UNIT_NAMES)
        specs.append({"shared": shared_df, "unit_specific": us})
    return pp.PanelParameters(theta=specs)


def _run_mif_batch_carry(pomp_dict, starts_batch, carry_batch, specific_names,
                         rw_sd_dict, nmif, Mp, Np, Np_rep, mif_key, pf_key, vmap_chunk_size):
    panel = pp.PanelPomp(
        Pomp_dict=pomp_dict,
        theta=make_panel_parameters_with_carry(starts_batch, carry_batch, specific_names),
    )
    panel.mif(J=Mp, M=nmif,
              rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
              key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size)
    panel.pfilter(J=Np, reps=Np_rep, key=pf_key,
                  chunk_size=vmap_chunk_size if vmap_chunk_size else 1)
    ll = np.asarray(panel.results_history[-1].logLiks.values)
    assert ll.ndim == 3, f"expected (Nstarts,U,reps), got {ll.shape}"
    unit_est, unit_se = logmeanexp_se(ll, axis=2)
    panel_ll = np.nansum(unit_est, axis=1)
    panel_se = np.sqrt(np.nansum(unit_se ** 2, axis=1))
    rows = []
    for idx, td in enumerate(panel.theta._to_list()):
        rows.append({
            "loglik": float(panel_ll[idx]),
            "se": float(panel_se[idx]),
            "coef": panel_theta_to_coef(td, specific_names),
        })
    return rows


def run_mif_round_carry(pomp_dict, starts_df, carry, specific_names, rw_sd_dict,
                        nmif, Mp, Np, Np_rep, key_seed, vmap_chunk_size):
    starts_df = starts_df.reset_index(drop=True)
    n = len(starts_df)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    rows = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts_df.iloc[lo:lo + BATCH_SIZE]
        carry_batch = carry[lo:lo + len(batch)]
        print(f"    batch {bi + 1}/{n_batches} (starts {lo}..{lo + len(batch) - 1})", flush=True)
        rows.extend(_run_mif_batch_carry(
            pomp_dict, batch, carry_batch, specific_names, rw_sd_dict, nmif,
            Mp, Np, Np_rep, jax.random.fold_in(mif_key, bi),
            jax.random.fold_in(pf_key, bi), vmap_chunk_size))
    return rows


# ============================================================================
# 7. Output in collect_lrt.R's schema (list(ll, se, coef, Np, Mp, Np_rep,
#    block, Nmif)) written as an .rds via pyreadr.
# ============================================================================
def write_result_rds(path: Path, ll: float, se: float, coef: dict[str, float],
                     Np: int, Mp: int, Np_rep: int, nmif_rounds: tuple[int, int]) -> None:
    """Write list(ll, se, coef, Np, Mp, Np_rep, block, Nmif) as an .rds.

    collect_lrt.R reads res$ll (and the Np/Mp/Np_rep/block guard).  pyreadr
    cannot write an arbitrary R named-list, so we write a one-row data.frame
    whose columns hold every scalar field plus the coef vector spread across
    "coef.<name>" columns and Nmif.round1 / Nmif.round2.  collect_lrt.R uses
    res$ll, res$se, res$Np, res$Mp, res$Np_rep, res$block via `$`, which on a
    data.frame selects the column -> works unchanged.

    *** collect_lrt.R adaptation (only if you point it at these files): ***
    collect_lrt.R does `res <- readRDS(f); res$ll`.  With this data.frame that
    returns the column (a length-1 vector) -> identical numeric value, so
    NO change is needed for $ll/$se/$Np/$Mp/$Np_rep/$block.  If you prefer a
    native R list, run the 3-line reshaper in the module docstring's
    VALIDATION section, or set PYPOMP_OVERWRITE=1 and re-save via R:
        x <- readRDS(f); saveRDS(as.list(x), f)   # then coef.* stay flat
    """
    import pyreadr
    path.parent.mkdir(parents=True, exist_ok=True)
    rec = {
        "ll": [float(ll)],
        "se": [float(se)],
        "Np": [int(Np)],
        "Mp": [int(Mp)],
        "Np_rep": [int(Np_rep)],
        "block": [True],
        "Nmif.round1": [int(nmif_rounds[0])],
        "Nmif.round2": [int(nmif_rounds[1])],
    }
    for k, v in coef.items():
        rec[f"coef.{k}"] = [float(v)]
    pyreadr.write_rds(str(path), pd.DataFrame(rec))


# ============================================================================
# 8. Main
# ============================================================================
def parse_args(argv):
    p = argparse.ArgumentParser(
        description="GPU/pypomp bootstrap-LRT fit for one (b, target) of SIRJPF2.")
    p.add_argument("b", type=int, help="bootstrap dataset index 1..100")
    p.add_argument("target", type=str,
                   help='"null" (all-shared) or an alternative name: '
                        + ", ".join(ALT_ORDER))
    p.add_argument("--run-level", type=int,
                   default=int(os.environ.get("PYPOMP_RUN_LEVEL", "3")),
                   choices=[1, 2, 3], help="particle/iteration budget (default 3 = R settings)")
    p.add_argument("--n-starts", type=int,
                   default=int(os.environ["PYPOMP_NSTARTS"]) if os.environ.get("PYPOMP_NSTARTS") else None,
                   help="mif starts (default = run-level; R uses 3*100=300 at level 3)")
    return p.parse_args(argv)


def main(argv) -> int:
    args = parse_args(argv)
    b = args.b
    target = args.target
    rl = args.run_level
    if not (1 <= b <= 100):
        raise ValueError("b must be 1..100")
    is_null = (target == "null")
    if not is_null and target not in ALT_PARAMS_MAP:
        raise ValueError(f"unknown target {target!r}; valid: null, " + ", ".join(ALT_ORDER))

    specific_names = [] if is_null else ALT_PARAMS_MAP[target]
    Mp = ALGO["Mp"][rl]
    Np = ALGO["Np"][rl]
    Np_rep = ALGO["Np_rep"][rl]
    nmif_rounds = ALGO["Nmif"][rl]
    nstarts = args.n_starts if args.n_starts is not None else ALGO["nstarts"][rl]
    vmap_chunk = os.environ.get("PYPOMP_VMAP_CHUNK")
    vmap_chunk = int(vmap_chunk) if vmap_chunk else None

    # seed offset: null vs alt index (stable, reproducible)  [R: set.seed(801000+b[+100000])]
    tgt_idx = 0 if is_null else (ALT_ORDER.index(target) + 1)

    print(f"=== pypomp bootstrap-LRT  b={b}  target={target} (run_level {rl}) ===")
    print(f"JAX backend = {jax.default_backend()}  x64 = {_USE_X64}  devices = {jax.devices()}")
    print(f"specific (unit-specific block) = {specific_names if specific_names else '(none; all-shared null)'}")
    print(f"Mp={Mp}  Np={Np}  Np_rep={Np_rep}  Nmif={nmif_rounds}  nstarts={nstarts}  batch={BATCH_SIZE}")

    base = Path(__file__).resolve().parent
    sim_data = read_sim_data(b)
    true_params = read_true_params()
    theta0 = {n: true_params[n] for n in PARAM_NAMES}
    pomp_dict = build_pomp_dict(sim_data, theta0)

    # ---- Round 1 starts: ALL warm starts at true_params  [R: shared.start=true_params,
    # specific.start=rep(true,8); replicated across the 3*cores foreach iterations].
    starts = pd.DataFrame([theta0] * nstarts)[PARAM_NAMES]

    print(f"Round 1 (Nmif={nmif_rounds[0]}, rw={RW_ROUNDS[0]})...", flush=True)
    round_one = run_mif_round(
        pomp_dict, starts, specific_names, make_rw_sd(RW_ROUNDS[0]),
        nmif=nmif_rounds[0], Mp=Mp, Np=Np, Np_rep=Np_rep,
        key_seed=801_000 + 1000 * tgt_idx + b, vmap_chunk_size=vmap_chunk)
    print(f"  best round-1 ll = {max(r['loglik'] for r in round_one):.3f}")

    # ---- Select top 25%, replicate x4 (carrying the unit-specific block) ----
    starts_df, carry = select_round_two(round_one, specific_names)

    print(f"Round 2 (Nmif={nmif_rounds[1]}, rw={RW_ROUNDS[1]})...", flush=True)
    round_two = run_mif_round_carry(
        pomp_dict, starts_df, carry, specific_names, make_rw_sd(RW_ROUNDS[1]),
        nmif=nmif_rounds[1], Mp=Mp, Np=Np, Np_rep=Np_rep,
        key_seed=801_000 + 100_000 + 1000 * tgt_idx + b, vmap_chunk_size=vmap_chunk)

    # ---- Best of round 2  [R: which.max(lls[1,]); ll, se, coef] ----
    best = max(range(len(round_two)), key=lambda i: round_two[i]["loglik"])
    ll = round_two[best]["loglik"]
    se = round_two[best]["se"]
    coef = round_two[best]["coef"]

    # ---- Output path  [R: results_null/lrt_null_<b>.rds | results_alt/lrt_<alt>_<b>.rds] ----
    overwrite = os.environ.get("PYPOMP_OVERWRITE", "0").lower() in {"1", "true", "yes"}
    if is_null:
        name = f"lrt_null_{b}" if overwrite else f"lrt_null_{b}_pypomp"
        out_path = base / "results_null" / f"{name}.rds"
    else:
        name = f"lrt_{target}_{b}" if overwrite else f"lrt_{target}_{b}_pypomp"
        out_path = base / "results_alt" / f"{name}.rds"

    write_result_rds(out_path, ll, se, coef, Np, Mp, Np_rep, nmif_rounds)

    print("\n================= RESULT =================")
    print(f"ll = {ll:.4f}   se = {se:.4f}   (compare to R {('results_null/lrt_null_%d.rds' % b) if is_null else ('results_alt/lrt_%s_%d.rds' % (target, b))})")
    print(f"coef entries = {len(coef)}  (shared {len([k for k in coef if '[' not in k])} + unit-specific {len([k for k in coef if '[' in k])})")
    print(f"saved -> {out_path}")
    print("==========================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
