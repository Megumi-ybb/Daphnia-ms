#!/usr/bin/env python3
"""
pypomp_coverage_profile.py
==========================
GPU / pypomp port of  coverage_profile.R  (the *general* version that takes the
profiled parameter as a command-line argument).  It reproduces the two-round
profile-MPIF coverage workflow used in the SIRJPF2 MCAP coverage study and writes
the SAME minimal output consumed by collect_coverage.R / calibrate_mcap.R:

    coverage_results/profile_<param>_<b>.rds

Usage
-----
    python pypomp_coverage_profile.py <dataset_index 1..100> <param>

    # examples (the four parameters of the coverage study + any profilable name)
    python pypomp_coverage_profile.py 42 rn
    python pypomp_coverage_profile.py 42 ri
    python pypomp_coverage_profile.py 42 sigF
    python pypomp_coverage_profile.py 42 sigP

Environment toggles (read BEFORE the heavy work; defaults reproduce the R run):
    PYPOMP_USE_GPU = 1 (default) | 0/cpu     -> JAX device
    PYPOMP_X64     = 1 (default) | 0          -> float64 (R uses double; keep 1)
    PYPOMP_RUN_LEVEL = 3 (default) | 1 | 2    -> particle/iteration budget
    PYPOMP_VMAP_CHUNK = (unset) | <int>       -> chunk the per-UNIT vmap inside pypomp
    PYPOMP_BATCH      = 400 (default) | <int>  -> outer batch over profile STARTS; this is
                                                 the knob that actually caps GPU memory,
                                                 because pypomp's vmap_chunk_size only chunks
                                                 the unit axis, NOT the 6400-wide start axis.
                                                 No effect on results (each start is an
                                                 independent fit, exactly as in R's foreach).

REQUIRES the development pypomp (>=0.4.6, with PanelPomp/PanelParameters/RWSigma/ParTrans;
`pip install -e` it) and a CUDA JAX (`pip install -U "jax[cuda12]"`).  Inputs are the
committed CSVs (sim_data_all.csv, true_params.csv) read with pandas, so NO R and NO rdata
are needed; the only I/O dependency is `pip install pyreadr` (to write the output .rds).

================================  FIDELITY NOTE  ===============================
JAX's PRNG is NOT the same generator as R's, so the particle draws and the mif
parameter perturbations are NOT bit-identical to a run of coverage_profile.R.
What this script reproduces EXACTLY is the *estimator and workflow*:

  * identical SIRJPF2 process / measurement model (term-for-term with the R
    Csnippets in coverage_profile.R lines 58-135),
  * the same profile design STRUCTURE -- the profiled parameter on a log grid over
    [true*0.1, true*10], crossed with nprof random log-uniform starts (nprof^2 rows),
    built in NumPy (NumPy RNG != R RNG, so the specific random starts differ),
  * the identical two-round MPIF structure, particle/iteration counts, cooling,
    per-parameter random-walk sds, top-25%-and-replicate selection, and the
    final group-by-profiled-parameter / keep-max-loglik output.

Consequently the two runs must be compared *statistically* (profile shape, MLE
location, MCAP CI agree within Monte-Carlo error) -- NOT byte-for-byte.
The R<->pypomp construct mapping is annotated inline as  [R: <line/desc>].
===============================================================================
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
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

# float64: R/pomp likelihoods are double precision.  States span 0..1e20 and the
# panel loglik sums 8 units x 10 obs, so single precision can drift materially.
# Keep x64 ON for fidelity (it is slower on consumer GPUs -- toggle with PYPOMP_X64=0).
_USE_X64 = os.environ.get("PYPOMP_X64", "1").lower() not in {"0", "false", "no"}
jax.config.update("jax_enable_x64", _USE_X64)

import jax.numpy as jnp  # noqa: E402
from jax.scipy.special import gammaln  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pypomp as pp  # noqa: E402


# ----------------------------------------------------------------------------
# 1. Constants  [R: coverage_profile.R lines 26-28, 46, 143-147, 222-227]
# ----------------------------------------------------------------------------
N_UNITS = 8
UNIT_NAMES = [f"u{i}" for i in range(1, N_UNITS + 1)]

# Run-level budgets, indexed by run_level in {1,2,3} (1-based, matching R).
# [R: algorithmic.params  lines 222-227];  Nmif is hard-coded to 200 (round 1)
# and 300 (round 2) in the R mif2() calls, NOT taken from algorithmic.params$Nmif.
ALGORITHMIC_PARAMS = {
    "Np": [50, 320, 1000],      # final pfilter particles      [R: $Np]
    "Np_rep": [2, 10, 20],      # pfilter replicates           [R: $Np_rep]
    "Mp": [50, 400, 500],       # mif2 particles               [R: $Mp]
}
NMIF_ROUND1 = 200               # [R: line 248  Nmif = 200]
NMIF_ROUND2 = 300               # [R: line 319  Nmif = 300]
COOLING_A = 0.7                 # [R: cooling.fraction.50 = 0.7]  (pypomp a == R cooling.fraction.50)

# Outer batch size over profile STARTS (the GPU-memory knob).  With nprof=80 the
# design has 6400 starts; pypomp vmaps the whole start axis at once, so we feed it
# in batches of BATCH_SIZE.  Each start is an independent fit (exactly as R's
# foreach over 6400 rows), so batching does NOT change results -- only peak memory.
BATCH_SIZE = int(os.environ.get("PYPOMP_BATCH", "400"))

# Profilable parameters in EXACTLY the R order; index gives the seed offset.
# [R: coverage_profile.R lines 26-28]
PROFILABLE = [
    "xi", "theta_Sn", "theta_Si", "theta_In", "theta_Ii", "theta_P", "theta_Jn", "theta_Ji",
    "rn", "ri", "f_Sn", "f_Si", "probn", "probi", "sigIn", "sigIi", "sigJn", "sigJi",
    "sigF", "sigP", "k_Sn", "k_Si", "k_In", "k_Ii",
]

# Canonical 26-parameter order  [R: param_names  lines 143-146]
PARAM_NAMES = [
    "xi", "sigSn", "sigIn", "sigSi", "sigIi", "sigF", "sigP",
    "theta_Sn", "theta_In", "theta_Si", "theta_P", "theta_Ii",
    "f_Sn", "f_Si", "rn", "ri", "probn", "probi",
    "k_Ii", "k_In", "k_Sn", "k_Si", "sigJi", "sigJn", "theta_Jn", "theta_Ji",
]

# State vector  [R: state_names  line 147];  order is irrelevant (dict-keyed).
SIRJPF_STATENAMES = [
    "Sn", "In", "Jn", "Si", "Ii", "Ji", "F", "P",
    "T_Sn", "T_In", "T_Si", "T_Ii", "error_count",
]

# Parameters carried on the LOG scale  [R: parameter_trans(log=...)  lines 137-141].
# NB sigSn, sigSi are fixed at 0 and are deliberately NOT log-transformed here
# (log(0) = -inf); they are pinned by rw_sd = 0, exactly as in R where rw.sd = 0
# keeps them at the boundary value 0.
SIRJPF_LOG_PARAMS = (
    "rn", "ri", "f_Sn", "f_Si", "probn", "probi", "xi",
    "theta_Sn", "theta_Si", "theta_In", "theta_Ii", "theta_Jn", "theta_Ji", "theta_P",
    "sigIn", "sigIi", "sigJn", "sigJi", "sigF", "sigP",
    "k_Sn", "k_Si", "k_In", "k_Ii",
)
FIXED_ZERO_PARAMS = ("sigSn", "sigSi")  # [R: sigSn = sigSi = 0, rw.sd = 0]


# ----------------------------------------------------------------------------
# 2. SIRJPF2 model in JAX  [R: Csnippets  coverage_profile.R lines 58-141]
#    (these match the R Csnippets term-for-term)
# ----------------------------------------------------------------------------
def sirjpf_rproc(X_, theta_, key, covars, t, dt):
    """One Euler-Maruyama step  [R: dyn_rpro, euler(delta.t = 1/4)]."""
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

    delta = 0.013          # [R: double delta = 0.013]
    mu_food = 0.37         # [R: + 0.37 * dt]
    lambda_J = 0.1         # [R: 0.1 * Jn ... maturation rate]

    keys = jax.random.split(key, 8)
    sqdt = jnp.sqrt(dt)
    noiSn = sigSn * sqdt * jax.random.normal(keys[0])
    noiIn = sigIn * sqdt * jax.random.normal(keys[1])
    noiSi = sigSi * sqdt * jax.random.normal(keys[2])
    noiIi = sigIi * sqdt * jax.random.normal(keys[3])
    noiJn = sigJn * sqdt * jax.random.normal(keys[4])
    noiJi = sigJi * sqdt * jax.random.normal(keys[5])
    noiF = sigF * sqdt * jax.random.normal(keys[6])
    noiP = sigP * sqdt * jax.random.normal(keys[7])

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
    # [R: F_term: note the native term uses (Sn + xi*In + 1*Jn), invasive (Si + xi*Ii + 1*Ji)]
    F_term = (-f_Sn * F * (Sn + xi * In + Jn) * dt
              - f_Si * F * (Si + xi * Ii + Ji) * dt
              - delta * F * dt + mu_food * dt + F * noiF)
    P_term = (30.0 * theta_In * In * dt + 30.0 * theta_Ii * Ii * dt
              - f_Sn * (Sn + xi * In) * P * dt - f_Si * (Si + xi * Ii) * P * dt
              - theta_P * P * dt - delta * P * dt + P * noiP)

    Sn_new = Sn + Sn_term
    In_new = In + In_term
    Jn_new = Jn + Jn_term
    Si_new = Si + Si_term
    Ii_new = Ii + Ii_term
    Ji_new = Ji + Ji_term
    F_new = F + F_term
    P_new = P + P_term
    # [R: if (t-4 ~ 0) P += 25]  -- parasite inoculation pulse at day 4
    P_new = P_new + jnp.where(jnp.abs(t - 4.0) < 0.001, 25.0, 0.0)

    # [R: soft constraints -> error_count, then states clamped to 0]
    def viol(x, hi):
        return (x < 0.0) | (x > hi)

    eps = 0.0
    eps += jnp.where(viol(Sn_new, 1e5), 1.0, 0.0)
    eps += jnp.where(viol(Si_new, 1e5), 1.0e6, 0.0)
    eps += jnp.where(viol(F_new, 1e20), 1.0e3, 0.0)
    eps += jnp.where(viol(In_new, 1e5), 1.0e-3, 0.0)
    eps += jnp.where(viol(Ii_new, 1e5), 1.0e-9, 0.0)
    eps += jnp.where(viol(Jn_new, 1e5), 1.0e-3, 0.0)
    eps += jnp.where(viol(Ji_new, 1e5), 1.0e-9, 0.0)
    eps += jnp.where(viol(P_new, 1e20) & (t > 3.9), 1.0e-6, 0.0)

    Sn_new = jnp.clip(Sn_new, 0.0, 1e5)
    In_new = jnp.clip(In_new, 0.0, 1e5)
    Jn_new = jnp.clip(Jn_new, 0.0, 1e5)
    Si_new = jnp.clip(Si_new, 0.0, 1e5)
    Ii_new = jnp.clip(Ii_new, 0.0, 1e5)
    Ji_new = jnp.clip(Ji_new, 0.0, 1e5)
    F_new = jnp.clip(F_new, 0.0, 1e20)
    P_new = jnp.clip(P_new, 0.0, 1e20)

    return {
        "Sn": Sn_new, "In": In_new, "Jn": Jn_new,
        "Si": Si_new, "Ii": Ii_new, "Ji": Ji_new,
        "F": F_new, "P": P_new,
        # [R: T_Sn = fabs(Sn); ... ]  measurement-ready non-negative densities
        "T_Sn": jnp.abs(Sn_new), "T_In": jnp.abs(In_new),
        "T_Si": jnp.abs(Si_new), "T_Ii": jnp.abs(Ii_new),
        # [R: error_count is an accumvar -> reset to 0 at each obs by pypomp]
        "error_count": error_count + eps,
    }


def sirjpf_rinit(theta_, key, covars, t0):
    """Deterministic initial state  [R: dyn_init  lines 107-116]."""
    return {
        "Sn": jnp.array(2.333), "In": jnp.array(0.0), "Jn": jnp.array(0.0),
        "Si": jnp.array(0.667), "Ii": jnp.array(0.0), "Ji": jnp.array(0.0),
        "F": jnp.array(16.667), "P": jnp.array(0.0),
        "T_Sn": jnp.array(0.0), "T_In": jnp.array(0.0),
        "T_Si": jnp.array(0.0), "T_Ii": jnp.array(0.0),
        "error_count": jnp.array(0.0),
    }


def _sirjpf_nb_logpmf(y, mu, size):
    """Negative-binomial log-pmf in (mu, size) form  [R: dnbinom_mu(y, size, mu)]."""
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    return (gammaln(y + size) - gammaln(size) - gammaln(y + 1.0)
            + size * jnp.log(size / (size + mu))
            + y * jnp.log(mu / (size + mu)))


def sirjpf_dmeas(Y_, X_, theta_, covars, t):
    """4-D NB measurement log-likelihood  [R: dmeas  lines 118-128]."""
    error_count = X_["error_count"]
    ll = (_sirjpf_nb_logpmf(Y_["dentadult"], X_["T_Sn"], theta_["k_Sn"])
          + _sirjpf_nb_logpmf(Y_["dentinf"], X_["T_In"], theta_["k_In"])
          + _sirjpf_nb_logpmf(Y_["lumadult"], X_["T_Si"], theta_["k_Si"])
          + _sirjpf_nb_logpmf(Y_["luminf"], X_["T_Ii"], theta_["k_Ii"]))
    # [R: if (error_count > 0.0) lik = -150]
    return jnp.where(error_count > 0.0, -150.0, ll)


def _sirjpf_nb_sample(key, mu, size):
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    k1, k2 = jax.random.split(key)
    g = jax.random.gamma(k1, size) * (mu / size)   # Gamma-Poisson = NB
    return jax.random.poisson(k2, g)


def sirjpf_rmeas(X_, theta_, key, covars, t):
    """Measurement simulator  [R: rmeas  lines 130-135]  (unused by profiling)."""
    keys = jax.random.split(key, 4)
    return jnp.array([
        _sirjpf_nb_sample(keys[0], X_["T_Sn"], theta_["k_Sn"]),
        _sirjpf_nb_sample(keys[1], X_["T_In"], theta_["k_In"]),
        _sirjpf_nb_sample(keys[2], X_["T_Si"], theta_["k_Si"]),
        _sirjpf_nb_sample(keys[3], X_["T_Ii"], theta_["k_Ii"]),
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
# 3. R-FREE I/O: read pre-converted CSVs with pandas (NO R, NO rdata), build the
#    profile design in NumPy, write the output .rds with `pyreadr`.
#    The CSVs (simulated_data/sim_data_all.csv, true_params.csv) are produced ONCE
#    from the original .rds and committed, so the GPU node needs no R-data reader.
#    Only dependency for I/O is `pyreadr` (output write):  pip install pyreadr
# ----------------------------------------------------------------------------
def read_sim_data(sim_dir: Path, b: int) -> dict[str, pd.DataFrame]:
    """Read the b-th simulated panel from sim_data_all.csv (cols: b,unit,day,4 obs)."""
    df = pd.read_csv(sim_dir / "sim_data_all.csv")
    df = df[df["b"] == b]
    if df.empty:
        raise ValueError(f"no rows for dataset b={b} in {sim_dir/'sim_data_all.csv'}")
    out = {}
    for unit_name in UNIT_NAMES:
        out[unit_name] = (df[df["unit"] == unit_name]
                          [["day", "dentadult", "dentinf", "lumadult", "luminf"]]
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


def generate_parameter_profile(true_params: dict[str, float], dataset_index: int,
                               prof_name: str, nprof: int = 80) -> pd.DataFrame:
    """
    Build the profile design in NumPy, reproducing pomp::profile_design's structure
    [R: coverage_profile.R lines 183-205]:
      * the profiled parameter takes `nprof` values on a log-uniform GRID over
        [true*0.1, true*10];
      * the other shared params (excluding sigSn, sigSi, profiled) get `nprof`
        random log-uniform starting points over the same per-parameter box;
      * the grid is crossed with the random starts -> nprof*nprof rows
        (6400 at nprof=80, matching the R design).  sigSn = sigSi = 0.
    NB: NumPy RNG != R RNG, so the *specific* random starts differ from a run of
    coverage_profile.R; the design is statistically equivalent (uniform multistart
    in the same box) -- consistent with the FIDELITY note.
    """
    pidx = PROFILABLE.index(prof_name) + 1               # 1-based, used in the seed
    seed = 805 * 1000 + 1000 * pidx + dataset_index      # reproducible (NumPy RNG)
    rng = np.random.default_rng(seed)
    lb = {n: true_params[n] * 0.1 for n in PARAM_NAMES}  # [R: shared_ub/100]
    ub = {n: true_params[n] * 10.0 for n in PARAM_NAMES}  # [R: true*10]
    grid = np.exp(np.linspace(np.log(lb[prof_name]), np.log(ub[prof_name]), nprof))
    others = [n for n in PARAM_NAMES if n not in ("sigSn", "sigSi", prof_name)]
    lo = np.array([np.log(lb[n]) for n in others])
    hi = np.array([np.log(ub[n]) for n in others])
    free_pts = np.exp(rng.uniform(lo, hi, size=(nprof, len(others))))   # nprof starts
    rows = []
    for g in grid:                                       # cross grid x free starts
        for r in range(nprof):
            row = {n: float(free_pts[r, j]) for j, n in enumerate(others)}
            row[prof_name] = float(g)
            row["sigSn"] = 0.0
            row["sigSi"] = 0.0
            rows.append(row)
    return pd.DataFrame(rows)[PARAM_NAMES]               # canonical column order


def write_profile_rds(profile_df: pd.DataFrame, path: Path) -> None:
    """Write the profile output as an .rds data.frame (read by collect_coverage.R)."""
    import pyreadr
    path.parent.mkdir(parents=True, exist_ok=True)
    pyreadr.write_rds(str(path), profile_df.reset_index(drop=True))


# ----------------------------------------------------------------------------
# 4. Per-parameter random-walk sds  [R: generate_sd  lines 207-218]
# ----------------------------------------------------------------------------
def generate_sd(value: float, profile_name: str) -> dict[str, float]:
    rw = {name: value for name in PARAM_NAMES}
    rw["sigSn"] = 0.0           # [R: sigSn = sigSi = 0]
    rw["sigSi"] = 0.0
    rw[profile_name] = 0.0       # [R: sd_list[profile_name] = 0]  (profiled param fixed)
    return rw


# ----------------------------------------------------------------------------
# 5. Panel construction + one MPIF round
# ----------------------------------------------------------------------------
def logmeanexp(values: np.ndarray, axis: int) -> np.ndarray:
    """log(mean(exp(x))) along `axis`, nan-robust."""
    m = np.nanmax(values, axis=axis, keepdims=True)
    out = np.squeeze(m, axis=axis) + np.log(np.nanmean(np.exp(values - m), axis=axis))
    return out


def build_pomp_dict(sim_data: dict[str, pd.DataFrame], true_params: dict[str, float]):
    """8 per-unit pp.Pomp objects  [R: pomplist loop  lines 153-175]."""
    theta = {name: true_params[name] for name in PARAM_NAMES}
    pomp_dict = {}
    for unit_name in UNIT_NAMES:
        ys_u = (sim_data[unit_name]
                .set_index("day")[["dentadult", "dentinf", "lumadult", "luminf"]]
                .astype(float))
        pomp_dict[unit_name] = pp.Pomp(
            ys=ys_u, theta=pp.PompParameters(theta), statenames=SIRJPF_STATENAMES, t0=1.0,
            rinit=sirjpf_rinit, rproc=sirjpf_rproc,
            dmeas=sirjpf_dmeas, rmeas=sirjpf_rmeas,
            par_trans=sirjpf_par_trans, dt=0.25,
            accumvars=("error_count",),       # [R: accumvars = c("error_count")]
        )
    return pomp_dict


def make_panel_parameters(param_rows: pd.DataFrame) -> "pp.PanelParameters":
    """
    Build an ALL-SHARED PanelParameters with one shared-theta spec per profile
    start  [R: panelPomp(pomplist, shared = ...)].  All 26 parameters are shared
    (including sigSn = sigSi = 0); unit_specific is empty.  Each row becomes one
    independent mif start -> the whole profile is vmapped on the GPU.
    """
    theta_specs = []
    empty_unit = pd.DataFrame(index=[], columns=UNIT_NAMES)
    for _, row in param_rows.iterrows():
        shared_df = pd.DataFrame(
            {"shared": [float(row[name]) for name in PARAM_NAMES]},
            index=PARAM_NAMES,
        )
        theta_specs.append({"shared": shared_df, "unit_specific": empty_unit})
    return pp.PanelParameters(theta=theta_specs)


def panel_theta_to_params(theta_dict: dict) -> dict[str, float]:
    """Extract the 26 shared coefficients from one panel-theta spec."""
    shared_df = theta_dict["shared"]
    params = {name: float(shared_df.loc[name, "shared"]) for name in shared_df.index}
    missing = [n for n in PARAM_NAMES if n not in params]
    if missing:
        raise ValueError(f"panel theta missing parameters: {missing}")
    return {name: params[name] for name in PARAM_NAMES}


def _run_mif_batch(pomp_dict, starts_batch, rw_sd_dict, nmif, mif_particles,
                   pfilter_particles, pfilter_reps, mif_key, pf_key, vmap_chunk_size):
    """MPIF + replicated pfilter on ONE batch of profile starts -> list of param dicts."""
    panel = pp.PanelPomp(Pomp_dict=pomp_dict, theta=make_panel_parameters(starts_batch))

    # [R: mif2(Nmif=nmif, rw.sd=..., cooling.fraction.50=0.7, Np=Mp)]
    # block=True is pypomp's MPIF; for an ALL-SHARED model it equals PIF (block is a
    # no-op without unit-specific params), so it matches the R all-shared mif2.
    # pypomp >=0.4.6: cooling is carried by RWSigma via .geometric_cooling(a=...),
    # not a mif(a=...) kwarg.  a=0.7 == R cooling.fraction.50=0.7.
    panel.mif(
        J=mif_particles, M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
        key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size,
    )

    # [R: ll <- replicate(Np_rep, unitLogLik(pfilter(m1, Np=Np)))]
    panel.pfilter(
        J=pfilter_particles, reps=pfilter_reps, key=pf_key,
        chunk_size=vmap_chunk_size if vmap_chunk_size else 1,
    )

    # logLiks: shape (Nstarts, U, reps) of per-unit, per-rep log-likelihoods.
    # [R: panel_logmeanexp(ll, MARGIN=1) = sum_units logmeanexp_over_reps].
    # IMPORTANT: logmeanexp is over reps PER UNIT, THEN summed over units --
    # NOT sum-then-logmeanexp (those do not commute).  Verified against pypomp's
    # (theta_idx, unit, rep) layout and panelPomp::panel_logmeanexp.
    loglik_array = np.asarray(panel.results_history[-1].logLiks.values)
    assert loglik_array.ndim == 3, f"expected (Nstarts,U,reps), got {loglik_array.shape}"
    per_unit_lme = logmeanexp(loglik_array, axis=2)      # (Nstarts, U)
    panel_loglik = np.nansum(per_unit_lme, axis=1)       # (Nstarts,)

    rows = []
    for idx, theta_dict in enumerate(panel.theta._to_list()):
        params = panel_theta_to_params(theta_dict)
        params["loglik"] = float(panel_loglik[idx])
        rows.append(params)
    return rows


def run_mif_round(pomp_dict, starts: pd.DataFrame, rw_sd_dict: dict[str, float],
                  nmif: int, mif_particles: int, pfilter_particles: int,
                  pfilter_reps: int, key_seed: int,
                  vmap_chunk_size) -> pd.DataFrame:
    """
    One round over ALL profile starts, fed to pypomp in batches of BATCH_SIZE.
    pypomp vmaps the entire start axis at once, so batching is the only way to bound
    the 6400-wide replicate dimension on a single GPU.  Each start is an independent
    fit (exactly like R's foreach over the 6400 rows), so the batched result equals
    an all-at-once result up to Monte-Carlo noise; batching only caps peak memory.
    Each batch gets an independent PRNG key via jax.random.fold_in.
    """
    starts = starts.reset_index(drop=True)
    n = len(starts)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    rows = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts.iloc[lo:lo + BATCH_SIZE]
        print(f"    batch {bi + 1}/{n_batches}  (starts {lo}..{lo + len(batch) - 1})", flush=True)
        rows.extend(_run_mif_batch(
            pomp_dict, batch, rw_sd_dict, nmif, mif_particles,
            pfilter_particles, pfilter_reps,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size,
        ))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# 6. Round-2 starts + final output  [R: lines 286-381]
# ----------------------------------------------------------------------------
def select_round_two_starts(round_one: pd.DataFrame) -> pd.DataFrame:
    """Top ceil(N/4) by loglik, each replicated 4x  [R: lines 290-298]."""
    n_select = math.ceil(len(round_one) / 4.0)
    selected = (round_one.sort_values("loglik", ascending=False)
                .head(n_select).drop(columns=["loglik"]).reset_index(drop=True))
    return selected.loc[selected.index.repeat(4)].reset_index(drop=True)


def make_output_profile(round_two: pd.DataFrame, profile_name: str) -> pd.DataFrame:
    """
    Keep, for each profiled-parameter value, the row with the MAX loglik, retaining
    ALL coefficient columns + loglik + log_<param>  [R: lines 362-381].
    (Full coefficients are required downstream for the composite ridge targets
    ri*f_Si and probi*f_Si in collect_coverage.R.)
    """
    idx = round_two.groupby(profile_name, sort=False)["loglik"].idxmax()
    best = round_two.loc[idx].reset_index(drop=True)
    best[f"log_{profile_name}"] = np.log(best[profile_name])
    return best


# ----------------------------------------------------------------------------
# 7. Main
# ----------------------------------------------------------------------------
def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="GPU/pypomp profile-coverage for one (parameter, dataset).")
    p.add_argument("dataset_index", type=int, help="dataset index 1..100")
    p.add_argument("param", type=str, help=f"profiled parameter; one of {PROFILABLE}")
    p.add_argument("--run-level", type=int,
                   default=int(os.environ.get("PYPOMP_RUN_LEVEL", "3")),
                   choices=[1, 2, 3], help="particle/iteration budget (default 3)")
    p.add_argument("--nprof", type=int, default=80, help="profile grid resolution (R: 80)")
    return p.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    b = args.dataset_index
    name_str = args.param
    rl = args.run_level
    if not (1 <= b <= 100):
        raise ValueError("dataset index must be 1..100")
    if name_str not in PROFILABLE:
        raise ValueError(f"unknown profile parameter {name_str!r}; valid: {PROFILABLE}")

    rl_i = rl - 1
    vmap_chunk = os.environ.get("PYPOMP_VMAP_CHUNK")
    vmap_chunk = int(vmap_chunk) if vmap_chunk else None

    print(f"=== pypomp coverage profile {name_str}: dataset {b} (run_level {rl}) ===")
    print(f"JAX backend = {jax.default_backend()}   x64 = {_USE_X64}   "
          f"devices = {jax.devices()}   vmap_chunk = {vmap_chunk}")

    base = Path(__file__).resolve().parent
    sim_dir = base / "simulated_data"
    out_path = base / "coverage_results" / f"profile_{name_str}_{b:03d}.rds"

    sim_data = read_sim_data(sim_dir, b)
    true_params = read_params(sim_dir)
    pomp_dict = build_pomp_dict(sim_data, true_params)
    parameter_shared = generate_parameter_profile(true_params, b, name_str, args.nprof)
    print(f"profile design rows = {len(parameter_shared)}")

    print("Round 1 starting...")
    round_one = run_mif_round(
        pomp_dict, parameter_shared, generate_sd(0.05, name_str),
        nmif=NMIF_ROUND1, mif_particles=ALGORITHMIC_PARAMS["Mp"][rl_i],
        pfilter_particles=ALGORITHMIC_PARAMS["Np"][rl_i],
        pfilter_reps=ALGORITHMIC_PARAMS["Np_rep"][rl_i],
        key_seed=10_000 + 1000 * (PROFILABLE.index(name_str) + 1) + b,
        vmap_chunk_size=vmap_chunk,
    )
    print("Round 1 done.")

    round_two_starts = select_round_two_starts(round_one)

    print("Round 2 starting...")
    round_two = run_mif_round(
        pomp_dict, round_two_starts, generate_sd(0.04, name_str),
        nmif=NMIF_ROUND2, mif_particles=ALGORITHMIC_PARAMS["Mp"][rl_i],
        pfilter_particles=ALGORITHMIC_PARAMS["Np"][rl_i],
        pfilter_reps=ALGORITHMIC_PARAMS["Np_rep"][rl_i],
        key_seed=20_000 + 1000 * (PROFILABLE.index(name_str) + 1) + b,
        vmap_chunk_size=vmap_chunk,
    )
    print("Round 2 done.")

    output_profile = make_output_profile(round_two, name_str)
    write_profile_rds(output_profile, out_path)

    print(f"Saved profile for {name_str} dataset {b}")
    print(f"  Rows: {len(output_profile)}   Columns: {list(output_profile.columns)}")
    print(f"  Path: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
