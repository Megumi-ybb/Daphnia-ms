#!/usr/bin/env python3
"""
pypomp_bootstrap_lrt.py  --  Single-species Lum SIRJPF parametric-bootstrap-LRT
================================================================================
GPU / pypomp port of  fit_null.R  +  fit_alt.R  (this directory) for the
SINGLE-SPECIES Lum SIRJPF family.  One bootstrap-LRT fit (one dataset b, one
target) becomes one batched MPIF run on the GPU, exactly as the SIRJPF2 coverage
study port (pypomp_coverage_profile.py) did for the profile workflow.

Usage
-----
    python pypomp_bootstrap_lrt.py <b 1..100> <target>

      target = "null"                      -> all-shared MPIF fit  (fit_null.R)
      target in {xi, theta_Si, theta_Ii,   -> alternative: that block is made
                 theta_P, probi, ri, f_Si}    UNIT-SPECIFIC, rest shared (fit_alt.R)

    # examples
    python pypomp_bootstrap_lrt.py 1 null
    python pypomp_bootstrap_lrt.py 1 theta_Ii
    python pypomp_bootstrap_lrt.py 42 f_Si

It reads simulated_data/sim_data_all.csv (pre-converted from the per-dataset
.rds; see CONVERTER below), fits via two-round MPIF, replicated-pfilters, and
writes the result in collect_lrt.R's schema:

    results_null/lrt_null_<b>.rds                 (target = null)
    results_alt/lrt_<target>_<b>.rds              (target = alternative)

Each .rds is an R list with $ll, $se, $coef, $Np, $Mp, $Np_rep, $block, $Nmif
(named c(round1=150, round2=250)) -- so collect_lrt.R runs UNCHANGED.

Environment toggles (defaults reproduce the R run):
    PYPOMP_USE_GPU = 1 (default) | 0/cpu      -> JAX device
    PYPOMP_X64     = 1 (default) | 0           -> float64 (R is double; keep 1)
    PYPOMP_RUN_LEVEL = 3 (default) | 1 | 2     -> particle/iteration budget
    PYPOMP_VMAP_CHUNK = (unset) | <int>        -> chunk the per-UNIT vmap in pypomp
    PYPOMP_BATCH      = 36  (default) | <int>  -> outer batch over the MPIF STARTS

REQUIRES the development pypomp (>=0.4.6 with PanelPomp/PanelParameters/RWSigma/
ParTrans) and a CUDA JAX.  R-FREE on the GPU node: inputs are the committed CSVs
read with pandas; output .rds is written with `pyreadr` (pip install pyreadr).

================================  CONVERTER  ===================================
The per-dataset simulated_data/sim_data_<b>.rds are converted ONCE (on a machine
with R) to a single flat simulated_data/sim_data_all.csv with columns
  b, unit, day, dentadult, dentinf
and true_params.rds to true_params.csv (name, value).  The exact R used was:

  B <- 100
  param_names <- c("sigSi","sigIi","sigF","sigP","f_Si","ri","k_Ii","k_Si",
                   "sigJi","probi","theta_Si","theta_Ii","theta_P","theta_Ji","xi")
  tp <- readRDS("simulated_data/true_params.rds")[param_names]
  write.csv(data.frame(name=names(tp), value=as.numeric(tp)),
            "simulated_data/true_params.csv", row.names=FALSE)
  out <- list(); k <- 1
  for (b in 1:B) {
    f <- sprintf("simulated_data/sim_data_%03d.rds", b); if (!file.exists(f)) next
    sd <- readRDS(f)
    for (i in 1:9) {
      u <- paste0("u", i); d <- sd[[u]]; colnames(d) <- c("day","dentadult","dentinf")
      out[[k]] <- data.frame(b=b, unit=u, day=d$day,
                             dentadult=d$dentadult, dentinf=d$dentinf); k <- k + 1
    }
  }
  write.csv(do.call(rbind, out), "simulated_data/sim_data_all.csv", row.names=FALSE)

The real .rds files are NOT modified.

================================  FIDELITY NOTE  ===============================
JAX's PRNG is not R's, so this is NOT bit-identical; it reproduces the same
estimator/workflow (model term-for-term, two-round MPIF, identical particle/
iteration counts, cooling, rw.sd, top-25%-and-replicate selection).  Compare the
two runs STATISTICALLY (ll within Monte-Carlo error), not byte-for-byte.

The R<->pypomp construct mapping is annotated inline as  [R: fit_*.R line/desc].
================================================================================

R Csnippet -> JAX TERM-BY-TERM MAPPING  (this family, single-species Lum SIRJPF)
-------------------------------------------------------------------------------
DROPPED vs SIRJPF2 (the 2-species mixed model): the entire n-species (dent)
compartment set -- states Sn,Jn,In,T_Sn,T_In and params rn,f_Sn,probn,theta_Sn,
theta_In,theta_Jn,sigSn,sigIn,sigJn,k_Sn,k_In -- because this family models only
the i-species (Lum).  KEPT: Si,Ii,Ji,F,P (+T_Si,T_Ii,error_count) and the
i-species params.  So 15 params / 9 units / 2 observables (vs SIRJPF2's 26/8/4).
NB the data OBSERVABLE COLUMNS are named "dentadult"/"dentinf" in the simulated
data (a labelling artifact of the shared simulate template) but they are the
i-species adult/infected counts, measured against T_Si (k_Si) and T_Ii (k_Ii)
-- exactly as fit_null.R's dmeas does.  See the per-line citations below.

rproc  [R: fit_null.R lines 21-54 / fit_alt.R 55-88]:
  noiSi = sigSi*sqrt(dt)*N(0,1)             [R L26]      <-> noiSi (keys[0])
  noiIi = sigIi*sqrt(dt)*N(0,1)             [R L27]      <-> noiIi (keys[1])
  noiJi = sigJi*sqrt(dt)*N(0,1)             [R L28]      <-> noiJi (keys[2])
  noiF  = sigF *sqrt(dt)*N(0,1)             [R L29]      <-> noiF  (keys[3])
  noiP  = sigP *sqrt(dt)*N(0,1)             [R L30]      <-> noiP  (keys[4])
  Si_term = 0.1*Ji*dt - theta_Si*Si*dt - probi*f_Si*Si*P*dt
            - delta*Si*dt + Si*noiSi        [R L32]      <-> Si_term
  Ji_term = ri*f_Si*F*Si*dt - 0.1*Ji*dt - theta_Ji*Ji*dt
            - delta*Ji*dt + Ji*noiJi        [R L33]      <-> Ji_term
  Ii_term = probi*f_Si*Si*P*dt - theta_Ii*Ii*dt
            - delta*Ii*dt + Ii*noiIi        [R L34]      <-> Ii_term
  F_term  = F*noiF - f_Si*F*(Si + xi*Ii + 1*Ji)*dt
            - delta*F*dt + 0.37*dt          [R L35]      <-> F_term
  P_term  = 30*theta_Ii*Ii*dt - f_Si*(Si + xi*Ii)*P*dt
            - theta_P*P*dt - delta*P*dt + P*noiP  [R L36] <-> P_term
  F+=,Si+=,Ii+=,Ji+=,P+=                    [R L38-42]   <-> *_new = * + *_term
  if(|t-4|<1e-3) P+=25                       [R L44]      <-> +jnp.where(|t-4|<1e-3,25,0)
  soft constraints (set offending state to 0.0, bump error_count):
    Si <0 or >1e5  -> Si=0, error_count+=1          [R L46]  <-> set Si_new=0
    F  <0 or >1e20 -> F =0, error_count+=1000        [R L47]  <-> set F_new=0
    Ii <0 or >1e5  -> Ii=0, error_count+=0.001       [R L48]  <-> set Ii_new=0
    Ji <0 or >1e5  -> Ji=0, error_count+=0.001       [R L49]  <-> set Ji_new=0
    (P<0 or >1e20) & t>3.9 -> P=0, error_count+=1e-6 [R L50]  <-> set P_new=0
    NOTE: R sets the offending state to 0.0 (NOT clip-to-upper-bound).  This is
    the family-correct choice and differs from the SIRJPF2 coverage port, which
    clipped.  It is observationally irrelevant at the violating obs (error_count
    >0 => lik=-150), but it changes the within-interval trajectory, so we match R.
  T_Si = fabs(Si); T_Ii = fabs(Ii)          [R L52-53]   <-> jnp.abs(Si_new/Ii_new)
  delta=0.013                               [R L24]
  accumvar error_count reset at each obs by pypomp (accumvars=("error_count",)).

rinit  [R: fit_null.R lines 56-65]:
  Si=3; F=16.667; Ji=0; T_Si=0; T_Ii=0; Ii=0; error_count=0; P=0
  (NB Si=3 here, not SIRJPF2's 0.667.)

dmeas  [R: fit_null.R lines 67-77]:
  if (error_count>0) lik=-150
  else lik = dnbinom_mu(dentadult, k_Si, T_Si, log) + dnbinom_mu(dentinf, k_Ii, T_Ii, log)
  -> R's dnbinom_mu(x, size, mu): dentadult ~ NB(size=k_Si, mu=T_Si),
                                  dentinf   ~ NB(size=k_Ii, mu=T_Ii).

partrans  [R: fit_null.R lines 84-88]:  log = all 15 params EXCEPT sigSi (sigSi
  is in the R log list but its true value is exactly 0 and rw.sd(sigSi)=0 pins it
  at the boundary; log(0)=-inf, so we treat sigSi as a fixed-zero param -- carried
  on the natural scale and never perturbed -- exactly mirroring how SIRJPF2 pins
  sigSn/sigSi at 0).

rw.sd  [R: fit_null.R lines 143-147]:  every param perturbed at dent_rw.sd EXCEPT
  sigSi=0.  k_Si and k_Ii ARE perturbed.  (Round1 sd=0.05, round2 sd=0.04.)

ALT unit-specific construction  [R: fit_alt.R create_parameters lines 37-48,
panelPomp(shared=..., specific=...) line 157]:  the alt block (one of xi,
theta_Si, theta_Ii, theta_P, probi, ri, f_Si) is moved out of `shared` and given
one column per unit (init = the shared/true value, replicated across 9 units);
everything else stays shared.  In pypomp this is a PanelParameters spec whose
`unit_specific` DataFrame has that block's row(s) with 9 unit columns, and whose
`shared` DataFrame holds the remaining params.  panel.mif(block=True) is MPIF;
block matters only because unit_specific is non-empty for the alt (it is a no-op
for the all-shared null, exactly as the R `block = TRUE` comment notes).
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
N_UNITS = 9                                                # [R: for i in 1:9]
UNIT_NAMES = [f"u{i}" for i in range(1, N_UNITS + 1)]

# Run-level budgets, indexed by run_level in {1,2,3} (1-based, matching R).
# Run level 3 reproduces fit_null.R / fit_alt.R EXACTLY: Np=Mp=1500, Np_rep=10.
# [R: fit_null.R lines 125-127]
ALGORITHMIC_PARAMS = {
    "Np": [50, 500, 1500],      # final pfilter particles      [R: Np = 1500]
    "Np_rep": [2, 10, 10],      # pfilter replicates           [R: Np_rep = 10]
    "Mp": [50, 500, 1500],      # mif2 particles               [R: Mp = 1500]
}
NMIF_ROUND1 = 150               # [R: Nmif = 150]
NMIF_ROUND2 = 250               # [R: Nmif = 250]
RW_ROUND1 = 0.05                # [R: dent_rw.sd = 0.05]  (round 1)
RW_ROUND2 = 0.04                # [R: dent_rw.sd = 0.04]  (round 2)
COOLING_A = 0.7                 # [R: cooling.fraction.50 = 0.7]

# Number of MPIF starts.  R uses 3 * getDoParWorkers() (= 3 * 100 = 300) but the
# fits are independent and the workflow keeps the top 25% then x4, so the start
# count is a robustness knob, not a definitional one.  Default 36 (the GPU port
# convention; raise with --n-starts to match R's 300 if desired).
N_STARTS_DEFAULT = int(os.environ.get("PYPOMP_NSTARTS", "36"))

# Outer batch over MPIF starts (GPU-memory knob; batching does not change results
# since each start is an independent fit -- exactly like R's foreach).
BATCH_SIZE = int(os.environ.get("PYPOMP_BATCH", "36"))

# Canonical 15-parameter order  [R: fit_null.R param_names lines 90-92]
PARAM_NAMES = [
    "sigSi", "sigIi", "sigF", "sigP", "f_Si",
    "ri", "k_Ii", "k_Si", "sigJi", "probi",
    "theta_Si", "theta_Ii", "theta_P", "theta_Ji", "xi",
]

# State vector  [R: fit_null.R state_names line 93];  order irrelevant (dict-keyed).
SIRJPF_STATENAMES = [
    "Si", "Ii", "Ji", "F", "P", "T_Si", "T_Ii", "error_count",
]

# Parameters carried on the LOG scale  [R: parameter_trans(log=...) lines 84-88].
# sigSi is in the R log list but is fixed at 0 (rw.sd=0); log(0)=-inf, so we pin
# it on the natural scale and never perturb it (same pattern as SIRJPF2's sigSn/Si).
SIRJPF_LOG_PARAMS = (
    "sigIi", "sigF", "sigP", "f_Si", "ri", "k_Ii", "k_Si", "sigJi", "probi",
    "theta_Si", "theta_Ii", "theta_P", "theta_Ji", "xi",
)
FIXED_ZERO_PARAMS = ("sigSi",)   # [R: sigSi = 0, rw.sd(sigSi) = 0]

# Alternatives this family tests  [R: fit_alt.R alt_params_map lines 18-26].
ALT_PARAMS_MAP = {
    "xi":       ["xi"],
    "theta_Si": ["theta_Si"],
    "theta_Ii": ["theta_Ii"],
    "theta_P":  ["theta_P"],
    "probi":    ["probi"],
    "ri":       ["ri"],
    "f_Si":     ["f_Si"],
}


# ----------------------------------------------------------------------------
# 2. Single-species Lum SIRJPF model in JAX
#    (term-for-term with the R Csnippets; mapping documented in the docstring)
# ----------------------------------------------------------------------------
def sirjpf_rproc(X_, theta_, key, covars, t, dt):
    """One Euler-Maruyama step  [R: euler(dyn_rpro, delta.t = 1/4)]."""
    Si, Ii, Ji = X_["Si"], X_["Ii"], X_["Ji"]
    F, P = X_["F"], X_["P"]
    error_count = X_["error_count"]

    sigSi, sigIi = theta_["sigSi"], theta_["sigIi"]
    sigJi, sigF, sigP = theta_["sigJi"], theta_["sigF"], theta_["sigP"]
    f_Si, ri, probi = theta_["f_Si"], theta_["ri"], theta_["probi"]
    theta_Si, theta_Ii = theta_["theta_Si"], theta_["theta_Ii"]
    theta_P, theta_Ji = theta_["theta_P"], theta_["theta_Ji"]
    xi = theta_["xi"]

    delta = 0.013          # [R L24: double delta = 0.013]

    keys = jax.random.split(key, 5)
    sqdt = jnp.sqrt(dt)
    noiSi = sigSi * sqdt * jax.random.normal(keys[0])   # [R L26]
    noiIi = sigIi * sqdt * jax.random.normal(keys[1])   # [R L27]
    noiJi = sigJi * sqdt * jax.random.normal(keys[2])   # [R L28]
    noiF = sigF * sqdt * jax.random.normal(keys[3])     # [R L29]
    noiP = sigP * sqdt * jax.random.normal(keys[4])     # [R L30]

    # [R L32]
    Si_term = (0.1 * Ji * dt - theta_Si * Si * dt
               - probi * f_Si * Si * P * dt - delta * Si * dt + Si * noiSi)
    # [R L33]
    Ji_term = (ri * f_Si * F * Si * dt - 0.1 * Ji * dt
               - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi)
    # [R L34]
    Ii_term = (probi * f_Si * Si * P * dt - theta_Ii * Ii * dt
               - delta * Ii * dt + Ii * noiIi)
    # [R L35]
    F_term = (F * noiF - f_Si * F * (Si + xi * Ii + 1.0 * Ji) * dt
              - delta * F * dt + 0.37 * dt)
    # [R L36]
    P_term = (30.0 * theta_Ii * Ii * dt - f_Si * (Si + xi * Ii) * P * dt
              - theta_P * P * dt - delta * P * dt + P * noiP)

    Si_new = Si + Si_term      # [R L39]
    Ii_new = Ii + Ii_term      # [R L40]
    Ji_new = Ji + Ji_term      # [R L41]
    F_new = F + F_term         # [R L38]
    P_new = P + P_term         # [R L42]
    # [R L44: if (|t-4| < 1e-3) P += 25]  -- parasite inoculation pulse at day 4
    P_new = P_new + jnp.where(jnp.abs(t - 4.0) < 0.001, 25.0, 0.0)

    # [R L46-50: soft constraints -> bump error_count, set offending state to 0.0]
    def viol(x, hi):
        return (x < 0.0) | (x > hi)

    bad_Si = viol(Si_new, 1e5)                      # [R L46]
    bad_F = viol(F_new, 1e20)                       # [R L47]
    bad_Ii = viol(Ii_new, 1e5)                      # [R L48]
    bad_Ji = viol(Ji_new, 1e5)                      # [R L49]
    bad_P = viol(P_new, 1e20) & (t > 3.9)           # [R L50]

    eps = 0.0
    eps += jnp.where(bad_Si, 1.0, 0.0)              # [R L46: error_count += 1]
    eps += jnp.where(bad_F, 1000.0, 0.0)            # [R L47: error_count += 1000]
    eps += jnp.where(bad_Ii, 0.001, 0.0)            # [R L48: error_count += 0.001]
    eps += jnp.where(bad_Ji, 0.001, 0.0)            # [R L49: error_count += 0.001]
    eps += jnp.where(bad_P, 0.000001, 0.0)          # [R L50: error_count += 1e-6]

    # R sets the offending state to 0.0 (NOT clip-to-hi).  Match R exactly.
    Si_new = jnp.where(bad_Si, 0.0, Si_new)
    F_new = jnp.where(bad_F, 0.0, F_new)
    Ii_new = jnp.where(bad_Ii, 0.0, Ii_new)
    Ji_new = jnp.where(bad_Ji, 0.0, Ji_new)
    P_new = jnp.where(bad_P, 0.0, P_new)

    return {
        "Si": Si_new, "Ii": Ii_new, "Ji": Ji_new, "F": F_new, "P": P_new,
        # [R L52-53: T_Si = fabs(Si); T_Ii = fabs(Ii)]
        "T_Si": jnp.abs(Si_new), "T_Ii": jnp.abs(Ii_new),
        # accumvar -> reset to 0 at each obs by pypomp
        "error_count": error_count + eps,
    }


def sirjpf_rinit(theta_, key, covars, t0):
    """Deterministic initial state  [R: dyn_init lines 56-65]."""
    return {
        "Si": jnp.array(3.0), "Ii": jnp.array(0.0), "Ji": jnp.array(0.0),
        "F": jnp.array(16.667), "P": jnp.array(0.0),
        "T_Si": jnp.array(0.0), "T_Ii": jnp.array(0.0),
        "error_count": jnp.array(0.0),
    }


def _nb_logpmf(y, mu, size):
    """NB log-pmf in (mu, size) form  [R: dnbinom_mu(y, size=size, mu=mu)]."""
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    return (gammaln(y + size) - gammaln(size) - gammaln(y + 1.0)
            + size * jnp.log(size / (size + mu))
            + y * jnp.log(mu / (size + mu)))


def sirjpf_dmeas(Y_, X_, theta_, covars, t):
    """2-D NB measurement log-likelihood  [R: dmeas lines 67-77].
    dentadult ~ NB(size=k_Si, mu=T_Si);  dentinf ~ NB(size=k_Ii, mu=T_Ii)."""
    ll = (_nb_logpmf(Y_["dentadult"], X_["T_Si"], theta_["k_Si"])
          + _nb_logpmf(Y_["dentinf"], X_["T_Ii"], theta_["k_Ii"]))
    # [R L68: if (error_count > 0.0) lik = -150]
    return jnp.where(X_["error_count"] > 0.0, -150.0, ll)


def _nb_sample(key, mu, size):
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    k1, k2 = jax.random.split(key)
    g = jax.random.gamma(k1, size) * (mu / size)    # Gamma-Poisson = NB
    return jax.random.poisson(k2, g)


def sirjpf_rmeas(X_, theta_, key, covars, t):
    """Measurement simulator  [R: rmeas lines 79-82]  (unused by the LRT fits)."""
    keys = jax.random.split(key, 2)
    return jnp.array([
        _nb_sample(keys[0], X_["T_Si"], theta_["k_Si"]),
        _nb_sample(keys[1], X_["T_Ii"], theta_["k_Ii"]),
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
# 3. R-FREE I/O: read the pre-converted CSVs with pandas; write .rds with pyreadr
# ----------------------------------------------------------------------------
def read_sim_data(sim_dir: Path, b: int) -> dict[str, pd.DataFrame]:
    """Read the b-th simulated panel from sim_data_all.csv
    (cols: b, unit, day, dentadult, dentinf)."""
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
    """Read true_params.csv (name, value)  [R: simulated_data/true_params.rds]."""
    tp = pd.read_csv(sim_dir / "true_params.csv")
    params = dict(zip(tp["name"].astype(str), tp["value"].astype(float)))
    missing = [n for n in PARAM_NAMES if n not in params]
    if missing:
        raise ValueError(f"true_params.csv missing parameters: {missing}")
    return {n: params[n] for n in PARAM_NAMES}


def write_result_rds(result: dict, path: Path) -> None:
    """Write the LRT result list as an .rds CONSUMED BY collect_lrt.R UNCHANGED.

    collect_lrt.R reads res$ll, res$se, res$coef, res$Np, res$Mp, res$Np_rep,
    res$block.  pyreadr writes a data.frame (not an R named list), so we write a
    single-row data.frame whose columns ARE those names: $ll, $se, $Np, $Mp,
    $Np_rep, $block, plus one column per coefficient prefixed `coef.<name>` and
    `Nmif.round1`/`Nmif.round2`.

    collect_lrt.R uses res$ll, res$se (scalars: a 1-col data.frame indexes the
    same as a scalar in R), res$Np/$Mp/$Np_rep/$block (scalars, used only for the
    paired-settings guard), and -- only when the obs fit is rebuilt -- res$coef.
    For the BOOTSTRAP arms collect_lrt.R touches ONLY $ll, $Np, $Mp, $Np_rep,
    $block, so this schema is sufficient and the guard sees matching settings.
    If a future caller needs the named-vector `coef`, re-save via the R one-liner
    in the header from these coef.<name> columns.
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
        "Nmif.round1": [int(result["Nmif"][0])],
        "Nmif.round2": [int(result["Nmif"][1])],
    }
    for name, val in result["coef"].items():
        row[f"coef.{name}"] = [float(val)]
    pyreadr.write_rds(str(path), pd.DataFrame(row))


# ----------------------------------------------------------------------------
# 4. rw.sd  [R: rw_sd(...) lines 143-147] -- ALL params perturbed except sigSi
# ----------------------------------------------------------------------------
def make_rw_sd(value: float) -> dict[str, float]:
    rw = {name: value for name in PARAM_NAMES}
    rw["sigSi"] = 0.0           # [R: sigSi = 0]
    return rw


# ----------------------------------------------------------------------------
# 5. Panel construction + one MPIF round
# ----------------------------------------------------------------------------
def logmeanexp_se(values: np.ndarray, axis: int):
    """logmeanexp and its Monte-Carlo SE along `axis`  [matches pomp::logmeanexp].
    [R: panel_logmeanexp(ll, MARGIN=1, se=TRUE)]."""
    m = np.nanmax(values, axis=axis, keepdims=True)
    w = np.exp(values - m)                       # weights in [0,1]
    n = np.sum(~np.isnan(values), axis=axis)
    mean_w = np.nanmean(w, axis=axis)
    est = np.squeeze(m, axis=axis) + np.log(mean_w)
    se = np.nanstd(w, axis=axis, ddof=1) / mean_w / np.sqrt(n)   # delta-method SE
    return est, se


def build_pomp_dict(sim_data: dict[str, pd.DataFrame], theta: dict[str, float]):
    """9 per-unit pp.Pomp objects  [R: pomplist loop lines 98-120]."""
    pomp_dict = {}
    for unit_name in UNIT_NAMES:
        ys_u = (sim_data[unit_name]
                .set_index("day")[["dentadult", "dentinf"]].astype(float))
        pomp_dict[unit_name] = pp.Pomp(
            ys=ys_u, theta=pp.PompParameters(theta), statenames=SIRJPF_STATENAMES,
            t0=1.0, rinit=sirjpf_rinit, rproc=sirjpf_rproc,
            dmeas=sirjpf_dmeas, rmeas=sirjpf_rmeas,
            par_trans=sirjpf_par_trans, dt=0.25,
            accumvars=("error_count",),       # [R: accumvars = c("error_count")]
        )
    return pomp_dict


def make_panel_parameters(param_rows: pd.DataFrame,
                          specific_names: list[str]) -> "pp.PanelParameters":
    """
    Build a PanelParameters with one spec per MPIF start.

    null  (specific_names == []):  ALL 15 params shared, unit_specific empty.
      [R: panelPomp(pomplist, shared = true_params)]
    alt   (specific_names != []):  the alt block goes UNIT-SPECIFIC (one column
      per unit, init = the row's value replicated across 9 units), the rest stay
      shared.  [R: create_parameters + panelPomp(shared=..., specific=...)]
    """
    shared_names = [n for n in PARAM_NAMES if n not in specific_names]
    specs = []
    for _, row in param_rows.iterrows():
        shared_df = pd.DataFrame(
            {"shared": [float(row[n]) for n in shared_names]}, index=shared_names)
        if specific_names:
            unit_df = pd.DataFrame(
                {u: [float(row[n]) for n in specific_names] for u in UNIT_NAMES},
                index=specific_names)
        else:
            unit_df = pd.DataFrame(index=[], columns=UNIT_NAMES)
        specs.append({"shared": shared_df, "unit_specific": unit_df})
    return pp.PanelParameters(theta=specs)


def panel_theta_to_coef(theta_dict: dict, specific_names: list[str]) -> dict[str, float]:
    """
    Extract the fitted coefficients in collect_lrt.R's `coef()` schema:
      * shared params by bare name,
      * unit-specific block params as "<param>[u1]".."<param>[u9]".
    [R: panelPomp::coef() naming, e.g. "xi[u1]".."xi[u9]" + shared rest].
    Matches the SIRJPF2 reference port (does NOT collapse the alt block).
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


def panel_theta_to_shared_row(theta_dict: dict) -> dict[str, float]:
    """Shared params only, for the round-1 -> round-2 warm-start carry."""
    sd = theta_dict["shared"]
    return {str(n): float(sd.loc[n, "shared"]) for n in sd.index}


def _run_mif_batch(pomp_dict, starts_batch, rw_sd_dict, nmif, Mp, Np, Np_rep,
                   specific_names, mif_key, pf_key, vmap_chunk_size):
    """MPIF + replicated pfilter on ONE batch of starts -> list of param dicts."""
    panel = pp.PanelPomp(
        Pomp_dict=pomp_dict,
        theta=make_panel_parameters(starts_batch, specific_names))

    # [R: mif2(Nmif, rw.sd, cooling.fraction.50=0.7, Np=Mp, block=TRUE)]
    # block=True is pypomp's MPIF.  For the all-shared null it is a no-op (no
    # unit_specific params) == R's "block = TRUE  # no-op for the all-shared null".
    # For the alt it block-resamples the unit-specific params == R's specific_model_block.
    # pypomp >=0.4.6: cooling lives on the RWSigma via .geometric_cooling(a=...),
    # a=0.7 == R cooling.fraction.50=0.7.
    panel.mif(
        J=Mp, M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
        key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size,
    )

    # [R: ll <- replicate(Np_rep, unitLogLik(pfilter(m1, Np=Np)))]
    panel.pfilter(
        J=Np, reps=Np_rep, key=pf_key,
        chunk_size=vmap_chunk_size if vmap_chunk_size else 1,
    )

    # logLiks: (Nstarts, U, reps).  panel_logmeanexp(MARGIN=1): logmeanexp over
    # reps PER UNIT, THEN sum over units (these do NOT commute) -- [R: line 154].
    loglik_array = np.asarray(panel.results_history[-1].logLiks.values)
    assert loglik_array.ndim == 3, f"expected (Nstarts,U,reps), got {loglik_array.shape}"
    per_unit_lme, per_unit_se = logmeanexp_se(loglik_array, axis=2)   # (Nstarts, U)
    panel_loglik = np.nansum(per_unit_lme, axis=1)                    # (Nstarts,)
    panel_se = np.sqrt(np.nansum(per_unit_se ** 2, axis=1))           # combine unit SEs

    rows = []
    for idx, theta_dict in enumerate(panel.theta._to_list()):
        rows.append({
            "loglik": float(panel_loglik[idx]),
            "se": float(panel_se[idx]),
            "coef": panel_theta_to_coef(theta_dict, specific_names),     # full per-unit coef (output schema)
            "shared": panel_theta_to_shared_row(theta_dict),            # for warm-start carry
            "theta": theta_dict,                                        # raw spec (carries unit-specific)
        })
    return rows


def run_mif_round(pomp_dict, starts: pd.DataFrame, rw_sd_dict, nmif, Mp, Np,
                  Np_rep, specific_names, key_seed, vmap_chunk_size) -> list:
    """One round over ALL starts, fed to pypomp in batches of BATCH_SIZE."""
    starts = starts.reset_index(drop=True)
    n = len(starts)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    rows = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts.iloc[lo:lo + BATCH_SIZE]
        print(f"    batch {bi + 1}/{n_batches}  (starts {lo}..{lo + len(batch) - 1})",
              flush=True)
        rows.extend(_run_mif_batch(
            pomp_dict, batch, rw_sd_dict, nmif, Mp, Np, Np_rep, specific_names,
            jax.random.fold_in(mif_key, bi), jax.random.fold_in(pf_key, bi),
            vmap_chunk_size))
    return rows


def select_round_two(round_one: list, specific_names: list[str]):
    """
    Top ceil(N/4) by loglik, each replicated x4, carrying BOTH the shared warm
    start AND (for the alt) the FULL per-unit unit-specific block from round 1.
    [R: fit_alt.R top-25% select + rep(each=4), inheriting @specific column-by-column].
    Matches the SIRJPF2 reference port -- does NOT collapse the alt block to a
    single value (an under-optimized alt biases Lambda DOWN and can MASK a real
    rejection).
    Returns (starts_df, carry list aligned to starts_df rows).
    """
    order = sorted(range(len(round_one)), key=lambda i: round_one[i]["loglik"], reverse=True)
    n_select = math.ceil(len(round_one) / 4.0)
    sel = order[:n_select]
    starts_rows, carry = [], []
    for i in sel:
        shared = round_one[i]["shared"]
        td = round_one[i]["theta"]
        for _ in range(4):                       # replicate x4
            starts_rows.append(shared)
            carry.append(td)
    starts_df = pd.DataFrame(starts_rows)
    # ensure every PARAM_NAME column exists; the alt block names are absent from
    # `shared`, so fill those from the carried unit-specific u1 value (the columns
    # are then overwritten per-unit by `carry` in make_panel_parameters_with_carry).
    for n in PARAM_NAMES:
        if n not in starts_df.columns:
            starts_df[n] = [float(carry[r]["unit_specific"].loc[n, UNIT_NAMES[0]])
                            for r in range(len(carry))]
    return starts_df[PARAM_NAMES], carry


def make_panel_parameters_with_carry(starts_df, carry, specific_names):
    """Like make_panel_parameters, but the alt's unit-specific block is taken
    from the carried round-1 per-unit matrices (NOT collapsed to one value),
    so round-2 resumes the per-unit estimates  [R: fit_alt.R specific.start]."""
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


def _run_mif_batch_carry(pomp_dict, starts_batch, carry_batch, rw_sd_dict, nmif,
                         Mp, Np, Np_rep, specific_names, mif_key, pf_key, vmap_chunk_size):
    """Same as _run_mif_batch but round-2 starts carry the full per-unit alt block."""
    panel = pp.PanelPomp(
        Pomp_dict=pomp_dict,
        theta=make_panel_parameters_with_carry(starts_batch, carry_batch, specific_names))
    panel.mif(
        J=Mp, M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]).geometric_cooling(a=COOLING_A),
        key=mif_key, block=True, vmap_chunk_size=vmap_chunk_size,
    )
    panel.pfilter(
        J=Np, reps=Np_rep, key=pf_key,
        chunk_size=vmap_chunk_size if vmap_chunk_size else 1,
    )
    loglik_array = np.asarray(panel.results_history[-1].logLiks.values)
    assert loglik_array.ndim == 3, f"expected (Nstarts,U,reps), got {loglik_array.shape}"
    per_unit_lme, per_unit_se = logmeanexp_se(loglik_array, axis=2)
    panel_loglik = np.nansum(per_unit_lme, axis=1)
    panel_se = np.sqrt(np.nansum(per_unit_se ** 2, axis=1))
    rows = []
    for idx, theta_dict in enumerate(panel.theta._to_list()):
        rows.append({
            "loglik": float(panel_loglik[idx]),
            "se": float(panel_se[idx]),
            "coef": panel_theta_to_coef(theta_dict, specific_names),
        })
    return rows


def run_mif_round_carry(pomp_dict, starts_df, carry, rw_sd_dict, nmif, Mp, Np,
                        Np_rep, specific_names, key_seed, vmap_chunk_size) -> list:
    """Round-2 driver: propagates the per-unit alt-block starts from round 1."""
    starts_df = starts_df.reset_index(drop=True)
    n = len(starts_df)
    mif_key = jax.random.key(key_seed)
    pf_key = jax.random.key(key_seed + 100_000)
    n_batches = math.ceil(n / BATCH_SIZE)
    rows = []
    for bi, lo in enumerate(range(0, n, BATCH_SIZE)):
        batch = starts_df.iloc[lo:lo + BATCH_SIZE]
        carry_batch = carry[lo:lo + len(batch)]
        print(f"    batch {bi + 1}/{n_batches}  (starts {lo}..{lo + len(batch) - 1})",
              flush=True)
        rows.extend(_run_mif_batch_carry(
            pomp_dict, batch, carry_batch, rw_sd_dict, nmif, Mp, Np, Np_rep,
            specific_names, jax.random.fold_in(mif_key, bi),
            jax.random.fold_in(pf_key, bi), vmap_chunk_size))
    return rows


# ----------------------------------------------------------------------------
# 6. Main: two-round MPIF  [R: fit_null.R / fit_alt.R]
# ----------------------------------------------------------------------------
def parse_args(argv: list[str]) -> argparse.Namespace:
    valid = ["null"] + list(ALT_PARAMS_MAP.keys())
    p = argparse.ArgumentParser(
        description="GPU/pypomp bootstrap-LRT for one (dataset, target), Lum SIRJPF.")
    p.add_argument("dataset_index", type=int, help="dataset index 1..100")
    p.add_argument("target", type=str,
                   help=f"'null' or an alternative; one of {valid}")
    p.add_argument("--run-level", type=int,
                   default=int(os.environ.get("PYPOMP_RUN_LEVEL", "3")),
                   choices=[1, 2, 3], help="particle/iteration budget (default 3 == R)")
    p.add_argument("--n-starts", type=int, default=N_STARTS_DEFAULT,
                   help="number of MPIF starts (R uses 3*100=300)")
    return p.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    b = args.dataset_index
    target = args.target
    rl = args.run_level
    if not (1 <= b <= 100):
        raise ValueError("dataset index must be 1..100")
    valid = ["null"] + list(ALT_PARAMS_MAP.keys())
    if target not in valid:
        raise ValueError(f"unknown target {target!r}; valid: {valid}")

    is_null = (target == "null")
    specific_names = [] if is_null else ALT_PARAMS_MAP[target]

    rl_i = rl - 1
    Mp = ALGORITHMIC_PARAMS["Mp"][rl_i]
    Np = ALGORITHMIC_PARAMS["Np"][rl_i]
    Np_rep = ALGORITHMIC_PARAMS["Np_rep"][rl_i]
    vmap_chunk = os.environ.get("PYPOMP_VMAP_CHUNK")
    vmap_chunk = int(vmap_chunk) if vmap_chunk else None

    print(f"=== pypomp Lum SIRJPF bootstrap-LRT: dataset {b}, target {target} "
          f"(run_level {rl}) ===")
    print(f"JAX backend = {jax.default_backend()}  x64 = {_USE_X64}  "
          f"devices = {jax.devices()}")
    print(f"Mp={Mp}  Np={Np}  Np_rep={Np_rep}  n_starts={args.n_starts}  "
          f"batch={BATCH_SIZE}  unit_specific={specific_names}")

    base = Path(__file__).resolve().parent
    sim_dir = base / "simulated_data"

    sim_data = read_sim_data(sim_dir, b)
    true_params = read_params(sim_dir)
    pomp_dict = build_pomp_dict(sim_data, true_params)

    # Round-1 starts: all = the warm start (true_params), exactly as fit_*.R
    # [R: shared.start = true_params (null) / parameter_candidates (alt)].
    theta0 = {n: true_params[n] for n in PARAM_NAMES}
    starts = pd.DataFrame([theta0] * args.n_starts)[PARAM_NAMES]

    # Seed offsets mirror the R set.seed scheme (null: 801*1000+b; alt: +100000).
    seed_base = 801_000 + b + (0 if is_null else 100_000)

    print("Round 1 starting...", flush=True)
    round_one = run_mif_round(
        pomp_dict, starts, make_rw_sd(RW_ROUND1),
        nmif=NMIF_ROUND1, Mp=Mp, Np=Np, Np_rep=Np_rep,
        specific_names=specific_names,
        key_seed=seed_base + 1000, vmap_chunk_size=vmap_chunk)
    print(f"Round 1 done.  best loglik = {max(r['loglik'] for r in round_one):.3f}",
          flush=True)

    # Select top 25%, replicate x4, carrying the FULL per-unit alt block forward.
    round_two_starts, carry = select_round_two(round_one, specific_names)

    print("Round 2 starting...", flush=True)
    round_two = run_mif_round_carry(
        pomp_dict, round_two_starts, carry, make_rw_sd(RW_ROUND2),
        nmif=NMIF_ROUND2, Mp=Mp, Np=Np, Np_rep=Np_rep,
        specific_names=specific_names,
        key_seed=seed_base + 2000, vmap_chunk_size=vmap_chunk)
    print(f"Round 2 done.  best loglik = {max(r['loglik'] for r in round_two):.3f}",
          flush=True)

    # [R: best <- which.max(lls[1,]); ll/se/coef from that mif]
    best = max(range(len(round_two)), key=lambda i: round_two[i]["loglik"])
    result = {
        "ll": float(round_two[best]["loglik"]),
        "se": float(round_two[best]["se"]),
        "coef": round_two[best]["coef"],     # per-unit coef ("<param>[u1]"..) for the alt
        "Np": Np, "Mp": Mp, "Np_rep": Np_rep, "block": True,
        "Nmif": (NMIF_ROUND1, NMIF_ROUND2),
    }

    # ---- Output path  [R: results_null/lrt_null_<b>.rds | results_alt/lrt_<alt>_<b>.rds] ----
    # PYPOMP_OVERWRITE=0 (default) writes a SIDECAR *_pypomp.rds (does NOT touch the
    # canonical R result); =1 overwrites the canonical lrt_*.rds for collect_lrt.R.
    overwrite = os.environ.get("PYPOMP_OVERWRITE", "0").lower() in {"1", "true", "yes"}
    if is_null:
        name = f"lrt_null_{b}" if overwrite else f"lrt_null_{b}_pypomp"
        out_path = base / "results_null" / f"{name}.rds"
    else:
        name = f"lrt_{target}_{b}" if overwrite else f"lrt_{target}_{b}_pypomp"
        out_path = base / "results_alt" / f"{name}.rds"

    write_result_rds(result, out_path)

    print("\n================= RESULT =================")
    print(f"ll = {result['ll']:.4f}   se = {result['se']:.4f}")
    coef = result["coef"]
    print(f"coef entries = {len(coef)}  (shared {len([k for k in coef if '[' not in k])} "
          f"+ unit-specific {len([k for k in coef if '[' in k])})")
    print(f"saved -> {out_path}")
    print("==========================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
