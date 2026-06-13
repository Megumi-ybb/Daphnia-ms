#!/usr/bin/env python3
"""
pypomp_all_shared.py
====================
GPU / pypomp replicate of  all_shared.R  (the canonical all-shared SIRJPF2 fit).

This is the VALIDATION GATE before the large GPU profile computation: it runs the
SAME model + the SAME marginalized panel iterated filtering (MPIF) on the SAME real
mesocosm data as all_shared.R, as a 3-round global search, and reports the best
log-likelihood and parameter estimate.  If this reproduces the R all-shared MLE
(published ll ~= -880.56, AIC 1809.12; mif.estimate ~= the `shared_parameter`
start), then pypomp/GPU is trustworthy and the profile (same model + algorithm,
just many more starts) can proceed.  See the memory note `project_pypomp_gpu_port`.

Usage
-----
    python pypomp_all_shared.py                 # run_level 3, 360 starts (= R 10*36 cores)
    python pypomp_all_shared.py --run-level 1   # quick smoke test
    python pypomp_all_shared.py --n-starts 40   # fewer starts for a fast check

Environment toggles (defaults reproduce the R run):
    PYPOMP_USE_GPU = 1 (default) | 0/cpu
    PYPOMP_X64     = 1 (default)               # R uses double precision
    PYPOMP_BATCH   = 120 (default) | <int>     # outer batch over starts (GPU memory)
    PYPOMP_VMAP_CHUNK = (unset) | <int>        # pypomp per-unit vmap chunk
    PYPOMP_DATA    = (auto) | <path to Mesocosmdata.xls>

REQUIRES the development pypomp tree (PanelPomp/PanelParameters/RWSigma/ParTrans),
NOT the old PyPI pypomp 0.0.4 -- put the dev pypomp on PYTHONPATH on the GPU host.
NO R is needed for this script: the .xls is read with pandas (`pip install xlrd` for
the legacy .xls format), and the result is written as CSV/JSON.

FIDELITY NOTE
-------------
JAX's PRNG is not R's, so this is NOT bit-identical; it reproduces the same
estimator/workflow.  The SIRJPF2 model functions below are byte-identical to
pypomp_coverage_profile.py (the profile port).  R<->pypomp mapping annotated as
[R: <all_shared.R line/desc>].
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

# --- JAX device / precision (must precede `import jax`) --------------------
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
# 1. Constants
# ----------------------------------------------------------------------------
N_UNITS = 8
UNIT_NAMES = [f"u{i}" for i in range(1, N_UNITS + 1)]   # u1=K, ..., u8=S  [R: trails]

# run-level budgets  [R: algorithmic.params lines 252-257];  Nmif is hard-coded to
# 150 (round 1), 150 (round 2), 400 (round 3) in the R mif2() calls.
ALGORITHMIC_PARAMS = {
    "Np": [50, 500, 1000],      # final pfilter particles
    "Np_rep": [2, 10, 10],      # pfilter replicates  (NOTE: 10, not 20 as in the profile)
    "Mp": [50, 500, 1000],      # mif2 particles
}
NMIF_ROUNDS = (150, 150, 400)   # [R: lines 277, 348, 420]
RW_ROUNDS = (0.05, 0.04, 0.04)  # [R: dent_rw.sd 0.05 then 0.04 then 0.04]
COOLING_A = 0.7                 # [R: cooling.fraction.50 = 0.7]
N_STARTS_DEFAULT = 360          # [R: 10 * getDoParWorkers() = 10 * 36]
BATCH_SIZE = int(os.environ.get("PYPOMP_BATCH", "120"))

# 26-parameter canonical order  [R: paramnames lines 228-229]
PARAM_NAMES = [
    "xi", "sigSn", "sigIn", "sigSi", "sigIi", "sigF", "sigP",
    "theta_Sn", "theta_In", "theta_Si", "theta_P", "theta_Ii",
    "f_Sn", "f_Si", "rn", "ri", "probn", "probi",
    "k_Ii", "k_In", "k_Sn", "k_Si", "sigJi", "sigJn", "theta_Jn", "theta_Ji",
]

# The all-shared MLE used both as panel coef and as the round-1 start
# [R: shared_parameter  lines 238-245]
SHARED_PARAMETER = {
    "ri": 13076.83, "rn": 59.04676, "f_Si": 1.838259e-05, "f_Sn": 0.001105668,
    "probi": 31.10083, "probn": 0.2565626, "xi": 28.6562,
    "theta_Sn": 0.1479834, "theta_Si": 0.0318604, "theta_Ii": 0.3531879,
    "theta_In": 0.5489315, "theta_P": 0.02024991, "theta_Ji": 0.0001299562,
    "theta_Jn": 0.0001532613, "sigSn": 0.0, "sigSi": 0.0, "sigIn": 0.0003063207,
    "sigIi": 0.02208698, "sigJi": 0.2727418, "sigJn": 0.2836891, "sigF": 0.1551729,
    "sigP": 0.238589, "k_Ii": 1.241092, "k_In": 1.005756, "k_Si": 4.715556,
    "k_Sn": 4.282648,
}

SIRJPF_STATENAMES = [
    "Sn", "In", "Jn", "Si", "Ii", "Ji", "F", "P",
    "T_Sn", "T_In", "T_Si", "T_Ii", "error_count",
]
SIRJPF_LOG_PARAMS = (
    "rn", "ri", "f_Sn", "f_Si", "probn", "probi", "xi",
    "theta_Sn", "theta_Si", "theta_In", "theta_Ii", "theta_Jn", "theta_Ji", "theta_P",
    "sigIn", "sigIi", "sigJn", "sigJi", "sigF", "sigP", "k_Sn", "k_Si", "k_In", "k_Ii",
)
FIXED_ZERO_PARAMS = ("sigSn", "sigSi")   # [R: sigSn = sigSi = 0, rw.sd = 0]


# ----------------------------------------------------------------------------
# 2. SIRJPF2 model in JAX  (identical to pypomp_coverage_profile.py and the
#    all_shared.R Csnippets, term-for-term)  [R: lines 41-204]
# ----------------------------------------------------------------------------
def sirjpf_rproc(X_, theta_, key, covars, t, dt):
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

    delta, mu_food, lambda_J = 0.013, 0.37, 0.1
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
    F_term = (-f_Sn * F * (Sn + xi * In + Jn) * dt
              - f_Si * F * (Si + xi * Ii + Ji) * dt
              - delta * F * dt + mu_food * dt + F * noiF)
    P_term = (30.0 * theta_In * In * dt + 30.0 * theta_Ii * Ii * dt
              - f_Sn * (Sn + xi * In) * P * dt - f_Si * (Si + xi * Ii) * P * dt
              - theta_P * P * dt - delta * P * dt + P * noiP)

    Sn_new = Sn + Sn_term; In_new = In + In_term; Jn_new = Jn + Jn_term
    Si_new = Si + Si_term; Ii_new = Ii + Ii_term; Ji_new = Ji + Ji_term
    F_new = F + F_term; P_new = P + P_term
    P_new = P_new + jnp.where(jnp.abs(t - 4.0) < 0.001, 25.0, 0.0)

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

    Sn_new = jnp.clip(Sn_new, 0.0, 1e5); In_new = jnp.clip(In_new, 0.0, 1e5)
    Jn_new = jnp.clip(Jn_new, 0.0, 1e5); Si_new = jnp.clip(Si_new, 0.0, 1e5)
    Ii_new = jnp.clip(Ii_new, 0.0, 1e5); Ji_new = jnp.clip(Ji_new, 0.0, 1e5)
    F_new = jnp.clip(F_new, 0.0, 1e20); P_new = jnp.clip(P_new, 0.0, 1e20)

    return {
        "Sn": Sn_new, "In": In_new, "Jn": Jn_new, "Si": Si_new, "Ii": Ii_new,
        "Ji": Ji_new, "F": F_new, "P": P_new,
        "T_Sn": jnp.abs(Sn_new), "T_In": jnp.abs(In_new),
        "T_Si": jnp.abs(Si_new), "T_Ii": jnp.abs(Ii_new),
        "error_count": error_count + eps,
    }


def sirjpf_rinit(theta_, key, covars, t0):
    return {
        "Sn": jnp.array(2.333), "In": jnp.array(0.0), "Jn": jnp.array(0.0),
        "Si": jnp.array(0.667), "Ii": jnp.array(0.0), "Ji": jnp.array(0.0),
        "F": jnp.array(16.667), "P": jnp.array(0.0),
        "T_Sn": jnp.array(0.0), "T_In": jnp.array(0.0),
        "T_Si": jnp.array(0.0), "T_Ii": jnp.array(0.0),
        "error_count": jnp.array(0.0),
    }


def _nb_logpmf(y, mu, size):
    mu = jnp.maximum(mu, 1e-10); size = jnp.maximum(size, 1e-10)
    return (gammaln(y + size) - gammaln(size) - gammaln(y + 1.0)
            + size * jnp.log(size / (size + mu)) + y * jnp.log(mu / (size + mu)))


def sirjpf_dmeas(Y_, X_, theta_, covars, t):
    ll = (_nb_logpmf(Y_["dentadult"], X_["T_Sn"], theta_["k_Sn"])
          + _nb_logpmf(Y_["dentinf"], X_["T_In"], theta_["k_In"])
          + _nb_logpmf(Y_["lumadult"], X_["T_Si"], theta_["k_Si"])
          + _nb_logpmf(Y_["luminf"], X_["T_Ii"], theta_["k_Ii"]))
    return jnp.where(X_["error_count"] > 0.0, -150.0, ll)


def _nb_sample(key, mu, size):
    mu = jnp.maximum(mu, 1e-10); size = jnp.maximum(size, 1e-10)
    k1, k2 = jax.random.split(key)
    return jax.random.poisson(k2, jax.random.gamma(k1, size) * (mu / size))


def sirjpf_rmeas(X_, theta_, key, covars, t):
    keys = jax.random.split(key, 4)
    return jnp.array([
        _nb_sample(keys[0], X_["T_Sn"], theta_["k_Sn"]),
        _nb_sample(keys[1], X_["T_In"], theta_["k_In"]),
        _nb_sample(keys[2], X_["T_Si"], theta_["k_Si"]),
        _nb_sample(keys[3], X_["T_Ii"], theta_["k_Ii"]),
    ], dtype=float)


def sirjpf_to_est(theta):
    out = {**theta}
    for n in SIRJPF_LOG_PARAMS:
        out[n] = jnp.log(jnp.maximum(theta[n], 1e-30))
    for n in FIXED_ZERO_PARAMS:
        out[n] = theta[n]
    return out


def sirjpf_from_est(theta):
    out = {**theta}
    for n in SIRJPF_LOG_PARAMS:
        out[n] = jnp.exp(theta[n])
    for n in FIXED_ZERO_PARAMS:
        out[n] = theta[n]
    return out


sirjpf_par_trans = pp.ParTrans(to_est=sirjpf_to_est, from_est=sirjpf_from_est)


# ----------------------------------------------------------------------------
# 3. Data: read Mesocosmdata.xls EXACTLY as all_shared.R does, via R
#    [R: lines 10-26]   (sheet 3, rows 91:170, day -> (day-1)*5+7, units K..S)
# ----------------------------------------------------------------------------
def _find_data_file() -> Path:
    if os.environ.get("PYPOMP_DATA"):
        return Path(os.environ["PYPOMP_DATA"])
    here = Path(__file__).resolve()
    for d in [here.parent, *here.parents]:
        for name in ("Mesocosmdata.xlsx", "Mesocosmdata.xls"):
            cand = d / name
            if cand.exists():
                return cand
    raise FileNotFoundError("Mesocosmdata.xls/.xlsx not found; set PYPOMP_DATA")


def read_mesocosm_data(data_path: Path) -> dict[str, pd.DataFrame]:
    """Read Mesocosmdata sheet 3 with PANDAS (no R), exactly as all_shared.R subsets it.
    [R: lines 10-26].  Reading a legacy .xls needs the `xlrd` engine: `pip install xlrd`
    (or convert Mesocosmdata.xls to .xlsx and point PYPOMP_DATA at it -> openpyxl)."""
    try:
        raw = pd.read_excel(data_path, sheet_name="both species combined")  # [R: sheet 3]
    except ImportError as e:
        raise ImportError(
            "Reading the legacy .xls needs xlrd:  pip install xlrd  "
            "(or convert Mesocosmdata.xls to .xlsx and set PYPOMP_DATA to it)."
        ) from e
    raw.columns = raw.columns.str.strip()
    cols = ["rep", "day", "dent.adult", "dent.inf", "lum.adult", "lum.adult.inf"]
    sub = raw.iloc[90:170][cols].copy()                 # [R: rows 91:170]
    sub["day"] = (sub["day"] - 1) * 5 + 7                # [R: day -> (day-1)*5+7]
    trails = ["K", "L", "M", "N", "O", "P", "Q", "S"]   # [R: trails -> u1..u8]
    obs = ["dent.adult", "dent.inf", "lum.adult", "lum.adult.inf"]
    out = {}
    for i, rep in enumerate(trails, start=1):
        d = sub[sub["rep"] == rep][["day"] + obs].copy()
        d.columns = ["day", "dentadult", "dentinf", "lumadult", "luminf"]
        out[f"u{i}"] = d.astype(float).sort_values("day")   # sort by day (== R reverse + order)
    return out


def save_result(mif_estimate: dict[str, float], loglik: float, se: float,
                base: Path) -> Path:
    """Save mif.estimate + loglik + se as CSV (+ JSON) -- no R required.
    Compare in R with:  read.csv('all_shared_pypomp.csv')  vs  all_shared.RData."""
    names = list(mif_estimate.keys()) + ["loglik", "se"]
    values = [mif_estimate[k] for k in mif_estimate] + [loglik, se]
    out_csv = base / "all_shared_pypomp.csv"
    pd.DataFrame({"name": names, "value": values}).to_csv(out_csv, index=False)
    (base / "all_shared_pypomp.json").write_text(
        json.dumps({"loglik": loglik, "se": se, "mif_estimate": mif_estimate}, indent=2))
    return out_csv


# ----------------------------------------------------------------------------
# 4. rw_sd  [R: rw_sd(...) lines 281-306]  -- ALL params perturbed except sigSn/sigSi
# ----------------------------------------------------------------------------
def make_rw_sd(value: float) -> dict[str, float]:
    rw = {name: value for name in PARAM_NAMES}
    rw["sigSn"] = 0.0
    rw["sigSi"] = 0.0
    return rw


# ----------------------------------------------------------------------------
# 5. Panel + one MPIF round (batched over starts for GPU memory)
# ----------------------------------------------------------------------------
def logmeanexp_se(values: np.ndarray, axis: int):
    """logmeanexp and its Monte-Carlo SE along `axis`  [matches pomp::logmeanexp]."""
    m = np.nanmax(values, axis=axis, keepdims=True)
    w = np.exp(values - m)                       # weights in [0,1]
    n = np.sum(~np.isnan(values), axis=axis)
    mean_w = np.nanmean(w, axis=axis)
    est = np.squeeze(m, axis=axis) + np.log(mean_w)
    se = np.nanstd(w, axis=axis, ddof=1) / mean_w / np.sqrt(n)   # delta-method SE
    return est, se


def build_pomp_dict(sim_data, theta):
    pomp_dict = {}
    for u in UNIT_NAMES:
        ys_u = (sim_data[u].set_index("day")[["dentadult", "dentinf", "lumadult", "luminf"]]
                .astype(float))
        pomp_dict[u] = pp.Pomp(
            ys=ys_u, theta=pp.PompParameters(theta), statenames=SIRJPF_STATENAMES, t0=1.0,
            rinit=sirjpf_rinit, rproc=sirjpf_rproc, dmeas=sirjpf_dmeas,
            rmeas=sirjpf_rmeas, par_trans=sirjpf_par_trans, dt=0.25,
            accumvars=("error_count",),
        )
    return pomp_dict


def make_panel_parameters(param_rows: pd.DataFrame):
    """All-shared PanelParameters, one shared spec per start (empty unit_specific)."""
    empty_unit = pd.DataFrame(index=[], columns=UNIT_NAMES)
    specs = []
    for _, row in param_rows.iterrows():
        shared_df = pd.DataFrame({"shared": [float(row[n]) for n in PARAM_NAMES]},
                                 index=PARAM_NAMES)
        specs.append({"shared": shared_df, "unit_specific": empty_unit})
    return pp.PanelParameters(theta=specs)


def panel_theta_to_params(theta_dict) -> dict[str, float]:
    sd = theta_dict["shared"]
    return {n: float(sd.loc[n, "shared"]) for n in PARAM_NAMES}


def _run_mif_batch(pomp_dict, starts_batch, rw_sd_dict, nmif, Mp, Np, Np_rep,
                   mif_key, pf_key, vmap_chunk_size):
    panel = pp.PanelPomp(Pomp_dict=pomp_dict, theta=make_panel_parameters(starts_batch))
    # [R: mif2(Nmif, rw.sd, cooling.fraction.50=0.7, Np=Mp)]  block=True == MPIF; all-shared => PIF.
    # pypomp >=0.4.6: cooling is carried by the RWSigma via .geometric_cooling(a=...),
    # NOT a mif(a=...) kwarg.  a=0.7 == R cooling.fraction.50=0.7.
    panel.mif(J=Mp, M=nmif,
              rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
              key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size)
    # [R: replicate(Np_rep, unitlogLik(pfilter(m1, Np))); panel_logmeanexp(MARGIN=1, se=TRUE)]
    panel.pfilter(J=Np, reps=Np_rep, key=pf_key,
                  chunk_size=vmap_chunk_size if vmap_chunk_size else 1)
    ll = np.asarray(panel.results_history[-1].logLiks.values)   # (Nstarts, U, reps)
    assert ll.ndim == 3, f"expected (Nstarts,U,reps), got {ll.shape}"
    unit_est, unit_se = logmeanexp_se(ll, axis=2)               # (Nstarts, U) each
    panel_ll = np.nansum(unit_est, axis=1)                      # sum over units
    panel_se = np.sqrt(np.nansum(unit_se ** 2, axis=1))         # combine unit SEs
    rows = []
    for idx, td in enumerate(panel.theta._to_list()):
        p = panel_theta_to_params(td)
        p["loglik"] = float(panel_ll[idx])
        p["se"] = float(panel_se[idx])
        rows.append(p)
    return rows


def run_mif_round(pomp_dict, starts, rw_sd_dict, nmif, Mp, Np, Np_rep,
                  key_seed, vmap_chunk_size) -> pd.DataFrame:
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
            pomp_dict, batch, rw_sd_dict, nmif, Mp, Np, Np_rep,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size))
    return pd.DataFrame(rows)


def select_and_replicate(df: pd.DataFrame) -> pd.DataFrame:
    """Top ceil(N/4) by loglik, drop loglik/se, replicate each 4x  [R: lines 327-332]."""
    n_select = math.ceil(len(df) / 4.0)
    sel = (df.sort_values("loglik", ascending=False).head(n_select)
           .drop(columns=["loglik", "se"]).reset_index(drop=True))
    return sel.loc[sel.index.repeat(4)].reset_index(drop=True)


# ----------------------------------------------------------------------------
# 6. Main: 3-round global search  [R: lines 266-468]
# ----------------------------------------------------------------------------
def parse_args(argv):
    p = argparse.ArgumentParser(description="GPU/pypomp replicate of all_shared.R")
    p.add_argument("--run-level", type=int, default=int(os.environ.get("PYPOMP_RUN_LEVEL", "3")),
                   choices=[1, 2, 3])
    p.add_argument("--n-starts", type=int, default=N_STARTS_DEFAULT,
                   help="number of search starts (R uses 10*36=360)")
    return p.parse_args(argv)


def main(argv) -> int:
    args = parse_args(argv)
    rl_i = args.run_level - 1
    Mp = ALGORITHMIC_PARAMS["Mp"][rl_i]
    Np = ALGORITHMIC_PARAMS["Np"][rl_i]
    Np_rep = ALGORITHMIC_PARAMS["Np_rep"][rl_i]
    vmap_chunk = os.environ.get("PYPOMP_VMAP_CHUNK")
    vmap_chunk = int(vmap_chunk) if vmap_chunk else None

    print(f"=== pypomp all_shared (run_level {args.run_level}, {args.n_starts} starts) ===")
    print(f"JAX backend = {jax.default_backend()}  x64 = {_USE_X64}  devices = {jax.devices()}")
    print(f"Mp(mif J) = {Mp}  Np(pfilter J) = {Np}  Np_rep = {Np_rep}  batch = {BATCH_SIZE}")

    base = Path(__file__).resolve().parent
    data_path = _find_data_file()
    print(f"data = {data_path}")
    sim_data = read_mesocosm_data(data_path)
    theta0 = {n: SHARED_PARAMETER[n] for n in PARAM_NAMES}
    pomp_dict = build_pomp_dict(sim_data, theta0)

    # Round 1 starts: all = the MLE  [R: parameter_candidates$shared, repeated]
    starts = pd.DataFrame([theta0] * args.n_starts)[PARAM_NAMES]

    for r in range(3):
        print(f"Round {r + 1} starting (Nmif={NMIF_ROUNDS[r]}, rw={RW_ROUNDS[r]})...")
        out = run_mif_round(pomp_dict, starts, make_rw_sd(RW_ROUNDS[r]),
                            nmif=NMIF_ROUNDS[r], Mp=Mp, Np=Np, Np_rep=Np_rep,
                            key_seed=805_000 + 1000 * (r + 1), vmap_chunk_size=vmap_chunk)
        best_ll = out["loglik"].max()
        print(f"Round {r + 1} done.  best loglik = {best_ll:.3f}")
        if r < 2:
            starts = select_and_replicate(out)

    # Best of round 3  [R: which.max(lls[1,]); mif.estimate; pf.loglik; se]
    best_row = out.loc[out["loglik"].idxmax()]
    mif_estimate = {n: float(best_row[n]) for n in PARAM_NAMES}
    pf_loglik = float(best_row["loglik"])
    pf_se = float(best_row["se"])

    out_csv = save_result(mif_estimate, pf_loglik, pf_se, base)

    print("\n================= RESULT (compare to all_shared.RData) =================")
    print(f"pf.loglik.of.mif.estimate = {pf_loglik:.4f}   (R all-shared published ll ~= -880.56)")
    print(f"s.e.                      = {pf_se:.4f}")
    print(f"saved mif.estimate + loglik + se -> {out_csv}")
    print("\n  param        start (MLE)        pypomp estimate")
    for n in PARAM_NAMES:
        print(f"  {n:<10} {SHARED_PARAMETER[n]:>16.6g}   {mif_estimate[n]:>16.6g}")
    print("========================================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
