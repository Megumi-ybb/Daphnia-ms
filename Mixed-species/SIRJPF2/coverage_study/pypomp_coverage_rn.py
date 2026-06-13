#!/usr/bin/env python3
"""
PyPOMP/GPU replicate of coverage_profile_rn.R.

Usage:
    python pypomp_coverage_rn.py <dataset_index>

This script reads the same simulated_data/sim_data_XXX.rds and
simulated_data/true_params.rds files as coverage_profile_rn.R, runs the same
two-round profile-MPIF workflow for rn, and writes the same minimal output:

    coverage_results/profile_rn_XXX.rds

The SIRJPF2 model functions are copied from the PyPOMP tutorial template in
PanelPomp-Python/daphnia_tut.qmd.
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

os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
if os.environ.get("PYPOMP_USE_GPU", "1").lower() in {"0", "false", "no", "cpu"}:
    os.environ["JAX_PLATFORMS"] = "cpu"
else:
    os.environ.pop("JAX_PLATFORMS", None)

import jax
import jax.numpy as jnp
from jax.scipy.special import gammaln
import numpy as np
import pandas as pd
import pypomp as pp


PROFILE_NAME = "rn"
RUN_LEVEL = 3
N_PROFILE = 80
N_UNITS = 8
UNIT_NAMES = [f"u{i}" for i in range(1, N_UNITS + 1)]

ALGORITHMIC_PARAMS = {
    "Np": [50, 320, 1000],
    "Np_rep": [2, 10, 20],
    "Mp": [50, 400, 500],
    "Nmif": [2, 320, 250],
}

PARAM_NAMES = [
    "xi",
    "sigSn",
    "sigIn",
    "sigSi",
    "sigIi",
    "sigF",
    "sigP",
    "theta_Sn",
    "theta_In",
    "theta_Si",
    "theta_P",
    "theta_Ii",
    "f_Sn",
    "f_Si",
    "rn",
    "ri",
    "probn",
    "probi",
    "k_Ii",
    "k_In",
    "k_Sn",
    "k_Si",
    "sigJi",
    "sigJn",
    "theta_Jn",
    "theta_Ji",
]

SIRJPF_STATENAMES = [
    "Sn",
    "In",
    "Jn",
    "Si",
    "Ii",
    "Ji",
    "F",
    "P",
    "T_Sn",
    "T_In",
    "T_Si",
    "T_Ii",
    "error_count",
]

SIRJPF_LOG_PARAMS = (
    "rn",
    "ri",
    "f_Sn",
    "f_Si",
    "probn",
    "probi",
    "xi",
    "theta_Sn",
    "theta_Si",
    "theta_In",
    "theta_Ii",
    "theta_Jn",
    "theta_Ji",
    "theta_P",
    "sigIn",
    "sigIi",
    "sigJn",
    "sigJi",
    "sigF",
    "sigP",
    "k_Sn",
    "k_Si",
    "k_In",
    "k_Ii",
)


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

    delta = 0.013
    mu_food = 0.37
    lambda_J = 0.1

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

    Sn_term = (
        lambda_J * Jn * dt
        - theta_Sn * Sn * dt
        - probn * f_Sn * Sn * P * dt
        - delta * Sn * dt
        + Sn * noiSn
    )
    Jn_term = (
        rn * f_Sn * F * Sn * dt
        - lambda_J * Jn * dt
        - theta_Jn * Jn * dt
        - delta * Jn * dt
        + Jn * noiJn
    )
    In_term = (
        probn * f_Sn * Sn * P * dt
        - theta_In * In * dt
        - delta * In * dt
        + In * noiIn
    )

    Si_term = (
        lambda_J * Ji * dt
        - theta_Si * Si * dt
        - probi * f_Si * Si * P * dt
        - delta * Si * dt
        + Si * noiSi
    )
    Ji_term = (
        ri * f_Si * F * Si * dt
        - lambda_J * Ji * dt
        - theta_Ji * Ji * dt
        - delta * Ji * dt
        + Ji * noiJi
    )
    Ii_term = (
        probi * f_Si * Si * P * dt
        - theta_Ii * Ii * dt
        - delta * Ii * dt
        + Ii * noiIi
    )

    F_term = (
        -f_Sn * F * (Sn + xi * In + Jn) * dt
        - f_Si * F * (Si + xi * Ii + Ji) * dt
        - delta * F * dt
        + mu_food * dt
        + F * noiF
    )
    P_term = (
        30.0 * theta_In * In * dt
        + 30.0 * theta_Ii * Ii * dt
        - f_Sn * (Sn + xi * In) * P * dt
        - f_Si * (Si + xi * Ii) * P * dt
        - theta_P * P * dt
        - delta * P * dt
        + P * noiP
    )

    Sn_new = Sn + Sn_term
    In_new = In + In_term
    Jn_new = Jn + Jn_term
    Si_new = Si + Si_term
    Ii_new = Ii + Ii_term
    Ji_new = Ji + Ji_term
    F_new = F + F_term
    P_new = P + P_term
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

    Sn_new = jnp.clip(Sn_new, 0.0, 1e5)
    In_new = jnp.clip(In_new, 0.0, 1e5)
    Jn_new = jnp.clip(Jn_new, 0.0, 1e5)
    Si_new = jnp.clip(Si_new, 0.0, 1e5)
    Ii_new = jnp.clip(Ii_new, 0.0, 1e5)
    Ji_new = jnp.clip(Ji_new, 0.0, 1e5)
    F_new = jnp.clip(F_new, 0.0, 1e20)
    P_new = jnp.clip(P_new, 0.0, 1e20)

    return {
        "Sn": Sn_new,
        "In": In_new,
        "Jn": Jn_new,
        "Si": Si_new,
        "Ii": Ii_new,
        "Ji": Ji_new,
        "F": F_new,
        "P": P_new,
        "T_Sn": jnp.abs(Sn_new),
        "T_In": jnp.abs(In_new),
        "T_Si": jnp.abs(Si_new),
        "T_Ii": jnp.abs(Ii_new),
        "error_count": error_count + eps,
    }


def sirjpf_rinit(theta_, key, covars, t0):
    return {
        "Sn": jnp.array(2.333),
        "In": jnp.array(0.0),
        "Jn": jnp.array(0.0),
        "Si": jnp.array(0.667),
        "Ii": jnp.array(0.0),
        "Ji": jnp.array(0.0),
        "F": jnp.array(16.667),
        "P": jnp.array(0.0),
        "T_Sn": jnp.array(0.0),
        "T_In": jnp.array(0.0),
        "T_Si": jnp.array(0.0),
        "T_Ii": jnp.array(0.0),
        "error_count": jnp.array(0.0),
    }


def _sirjpf_nb_logpmf(y, mu, size):
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    return (
        gammaln(y + size)
        - gammaln(size)
        - gammaln(y + 1.0)
        + size * jnp.log(size / (size + mu))
        + y * jnp.log(mu / (size + mu))
    )


def sirjpf_dmeas(Y_, X_, theta_, covars, t):
    error_count = X_["error_count"]
    ll_dent_adult = _sirjpf_nb_logpmf(Y_["dentadult"], X_["T_Sn"], theta_["k_Sn"])
    ll_dent_inf = _sirjpf_nb_logpmf(Y_["dentinf"], X_["T_In"], theta_["k_In"])
    ll_lum_adult = _sirjpf_nb_logpmf(Y_["lumadult"], X_["T_Si"], theta_["k_Si"])
    ll_lum_inf = _sirjpf_nb_logpmf(Y_["luminf"], X_["T_Ii"], theta_["k_Ii"])
    ll = ll_dent_adult + ll_dent_inf + ll_lum_adult + ll_lum_inf
    return jnp.where(error_count > 0.0, -150.0, ll)


def _sirjpf_nb_sample(key, mu, size):
    mu = jnp.maximum(mu, 1e-10)
    size = jnp.maximum(size, 1e-10)
    k1, k2 = jax.random.split(key)
    scale = mu / size
    g = jax.random.gamma(k1, size) * scale
    return jax.random.poisson(k2, g)


def sirjpf_rmeas(X_, theta_, key, covars, t):
    keys = jax.random.split(key, 4)
    y_dent_adult = _sirjpf_nb_sample(keys[0], X_["T_Sn"], theta_["k_Sn"])
    y_dent_inf = _sirjpf_nb_sample(keys[1], X_["T_In"], theta_["k_In"])
    y_lum_adult = _sirjpf_nb_sample(keys[2], X_["T_Si"], theta_["k_Si"])
    y_lum_inf = _sirjpf_nb_sample(keys[3], X_["T_Ii"], theta_["k_Ii"])
    return jnp.array([y_dent_adult, y_dent_inf, y_lum_adult, y_lum_inf], dtype=float)


def sirjpf_to_est(theta):
    out = {**theta}
    for name in SIRJPF_LOG_PARAMS:
        out[name] = jnp.log(jnp.maximum(theta[name], 1e-30))
    out["sigSn"] = theta["sigSn"]
    out["sigSi"] = theta["sigSi"]
    return out


def sirjpf_from_est(theta):
    out = {**theta}
    for name in SIRJPF_LOG_PARAMS:
        out[name] = jnp.exp(theta[name])
    out["sigSn"] = theta["sigSn"]
    out["sigSi"] = theta["sigSi"]
    return out


sirjpf_par_trans = pp.ParTrans(to_est=sirjpf_to_est, from_est=sirjpf_from_est)


def run_rscript(code: str, args: list[str]) -> str:
    command = ["Rscript", "--vanilla", "-e", code, *args]
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return completed.stdout


def read_sim_data_rds(path: Path) -> dict[str, pd.DataFrame]:
    code = r"""
    args <- commandArgs(trailingOnly = TRUE)
    x <- readRDS(args[1])
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("R package jsonlite is required to read RDS inputs from Python")
    }
    cat(jsonlite::toJSON(x, dataframe = "rows", auto_unbox = TRUE, digits = 17))
    """
    raw = run_rscript(code, [str(path)])
    parsed = json.loads(raw)
    out = {}
    for unit_name in UNIT_NAMES:
        unit_rows = parsed[unit_name]
        unit_df = pd.DataFrame(unit_rows)
        unit_df.columns = ["day", "dentadult", "dentinf", "lumadult", "luminf"]
        out[unit_name] = unit_df.astype(float).sort_values("day")
    return out


def read_params_rds(path: Path) -> dict[str, float]:
    code = r"""
    args <- commandArgs(trailingOnly = TRUE)
    x <- readRDS(args[1])
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("R package jsonlite is required to read RDS inputs from Python")
    }
    df <- data.frame(name = names(x), value = as.numeric(x), check.names = FALSE)
    cat(jsonlite::toJSON(df, dataframe = "rows", auto_unbox = TRUE, digits = 17))
    """
    raw = run_rscript(code, [str(path)])
    rows = json.loads(raw)
    params = {row["name"]: float(row["value"]) for row in rows}
    missing = [name for name in PARAM_NAMES if name not in params]
    if missing:
        raise ValueError(f"true_params.rds is missing parameters: {missing}")
    return params


def write_profile_rds(profile_df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        profile_df.to_csv(tmp_path, index=False)
        code = r"""
        args <- commandArgs(trailingOnly = TRUE)
        input_csv <- args[1]
        output_rds <- args[2]
        df <- read.csv(input_csv, check.names = FALSE)
        if (requireNamespace("tibble", quietly = TRUE)) {
          df <- tibble::as_tibble(df)
        }
        saveRDS(df, file = output_rds)
        """
        run_rscript(code, [str(tmp_path), str(path)])
    finally:
        tmp_path.unlink(missing_ok=True)


def generate_parameter_profile(params_path: Path, dataset_index: int) -> pd.DataFrame:
    code = r"""
    args <- commandArgs(trailingOnly = TRUE)
    true_params <- readRDS(args[1])
    dataset_index <- as.integer(args[2])
    prof_name <- args[3]
    nprof <- as.integer(args[4])
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("R package jsonlite is required to return profile starts to Python")
    }
    if (!requireNamespace("pomp", quietly = TRUE)) {
      stop("R package pomp is required to reproduce pomp::profile_design")
    }
    set.seed(802 * 1000 + dataset_index)
    shared_ub <- true_params * 10
    shared_lb <- shared_ub / 100
    ub_unit <- log(shared_ub[prof_name])
    lb_unit <- log(shared_lb[prof_name])
    shared_lb <- shared_lb[!(names(shared_lb) %in% c("sigSn", "sigSi", prof_name))]
    shared_ub <- shared_ub[!(names(shared_ub) %in% c("sigSn", "sigSi", prof_name))]
    parameter_shared <- pomp::profile_design(
      temp = seq(lb_unit, ub_unit, length.out = nprof),
      lower = log(shared_lb),
      upper = log(shared_ub),
      type = "runif",
      nprof = nprof
    )
    names(parameter_shared)[names(parameter_shared) == "temp"] <- prof_name
    parameter_shared <- exp(parameter_shared)
    parameter_shared$sigSn <- 0
    parameter_shared$sigSi <- 0
    cat(jsonlite::toJSON(parameter_shared, dataframe = "rows", auto_unbox = TRUE, digits = 17))
    """
    raw = run_rscript(code, [str(params_path), str(dataset_index), PROFILE_NAME, str(N_PROFILE)])
    rows = json.loads(raw)
    profile_df = pd.DataFrame(rows)
    missing = [name for name in PARAM_NAMES if name not in profile_df.columns]
    if missing:
        raise ValueError(f"profile design is missing parameters: {missing}")
    return profile_df


def generate_sd(value: float, profile_name: str) -> dict[str, float]:
    rw_sd = {
        "ri": value,
        "rn": value,
        "f_Si": value,
        "f_Sn": value,
        "probi": value,
        "probn": value,
        "xi": value,
        "theta_Sn": value,
        "theta_Si": value,
        "theta_Ii": value,
        "theta_In": value,
        "theta_P": value,
        "theta_Ji": value,
        "theta_Jn": value,
        "sigSn": 0.0,
        "sigSi": 0.0,
        "sigIn": value,
        "sigIi": value,
        "sigJi": value,
        "sigJn": value,
        "sigF": value,
        "sigP": value,
        "k_Ii": value,
        "k_In": value,
        "k_Si": value,
        "k_Sn": value,
    }
    rw_sd[profile_name] = 0.0
    return rw_sd


def build_pomp_dict(sim_data: dict[str, pd.DataFrame], true_params: dict[str, float]):
    theta = {name: true_params[name] for name in PARAM_NAMES}
    pomp_dict = {}
    for unit_name in UNIT_NAMES:
        ys_u = (
            sim_data[unit_name]
            .set_index("day")[["dentadult", "dentinf", "lumadult", "luminf"]]
            .astype(float)
        )
        pomp_dict[unit_name] = pp.Pomp(
            ys=ys_u,
            theta=theta,
            statenames=SIRJPF_STATENAMES,
            t0=1.0,
            rinit=sirjpf_rinit,
            rproc=sirjpf_rproc,
            dmeas=sirjpf_dmeas,
            rmeas=sirjpf_rmeas,
            par_trans=sirjpf_par_trans,
            dt=0.25,
            accumvars=("error_count",),
        )
    return pomp_dict


def make_panel_parameters(param_rows: pd.DataFrame) -> pp.PanelParameters:
    theta_specs = []
    shared_keys = [name for name in PARAM_NAMES if name not in {"sigSn", "sigSi"}]
    for _, row in param_rows.iterrows():
        shared_df = pd.DataFrame(
            {"shared": [float(row[name]) for name in shared_keys]},
            index=shared_keys,
        )
        unit_df = pd.DataFrame(
            {unit_name: [float(row["sigSn"]), float(row["sigSi"])] for unit_name in UNIT_NAMES},
            index=["sigSn", "sigSi"],
        )
        theta_specs.append({"shared": shared_df, "unit_specific": unit_df})
    return pp.PanelParameters(theta=theta_specs)


def panel_theta_to_params(theta_dict: dict[str, pd.DataFrame]) -> dict[str, float]:
    params = {}
    shared_df = theta_dict["shared"]
    for name in shared_df.index:
        params[name] = float(shared_df.loc[name, "shared"])

    unit_df = theta_dict["unit_specific"]
    for name in unit_df.index:
        values = unit_df.loc[name].dropna().astype(float)
        params[name] = float(values.iloc[0]) if len(values) else 0.0

    missing = [name for name in PARAM_NAMES if name not in params]
    if missing:
        raise ValueError(f"Panel theta is missing parameters: {missing}")
    return {name: params[name] for name in PARAM_NAMES}


def logmeanexp(values: np.ndarray, axis: int) -> np.ndarray:
    max_values = np.nanmax(values, axis=axis, keepdims=True)
    centered = np.exp(values - max_values)
    return np.squeeze(max_values, axis=axis) + np.log(np.nanmean(centered, axis=axis))


def run_mif_round(
    pomp_dict,
    starts: pd.DataFrame,
    rw_sd_dict: dict[str, float],
    nmif: int,
    mif_particles: int,
    pfilter_particles: int,
    pfilter_reps: int,
    key_seed: int,
) -> pd.DataFrame:
    panel = pp.PanelPomp(
        Pomp_dict=pomp_dict,
        theta=make_panel_parameters(starts),
    )
    panel.mif(
        J=mif_particles,
        M=nmif,
        rw_sd=pp.RWSigma(sigmas=rw_sd_dict, init_names=[]),
        a=0.7,
        key=jax.random.key(key_seed),
    )
    panel.pfilter(
        J=pfilter_particles,
        reps=pfilter_reps,
        key=jax.random.key(key_seed + 100_000),
    )

    loglik_array = np.asarray(panel.results_history[-1].logLiks.values)
    panel_loglik_by_rep = np.nansum(loglik_array, axis=1)
    panel_loglik = logmeanexp(panel_loglik_by_rep, axis=1)
    theta_list = panel.theta.to_list()

    rows = []
    for idx, theta_dict in enumerate(theta_list):
        params = panel_theta_to_params(theta_dict)
        params["loglik"] = float(panel_loglik[idx])
        rows.append(params)
    return pd.DataFrame(rows)


def select_round_two_starts(round_one: pd.DataFrame) -> pd.DataFrame:
    n_select = math.ceil(len(round_one) / 4.0)
    selected = round_one.sort_values("loglik", ascending=False).head(n_select)
    selected = selected.drop(columns=["loglik"]).reset_index(drop=True)
    return selected.loc[selected.index.repeat(4)].reset_index(drop=True)


def make_output_profile(round_two: pd.DataFrame) -> pd.DataFrame:
    best = (
        round_two.sort_values("loglik", ascending=True)
        .groupby(PROFILE_NAME, as_index=False)
        .tail(1)[[PROFILE_NAME, "loglik"]]
        .reset_index(drop=True)
    )
    best[f"log_{PROFILE_NAME}"] = np.log(best[PROFILE_NAME])
    return best


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run PyPOMP/GPU profile coverage for rn on one simulated dataset."
    )
    parser.add_argument("dataset_index", type=int, help="Dataset index between 1 and 100")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    dataset_index = args.dataset_index
    if dataset_index < 1 or dataset_index > 100:
        raise ValueError("Dataset index must be an integer between 1 and 100.")

    print(f"=== PyPOMP coverage profile {PROFILE_NAME}: dataset {dataset_index} ===")
    print(f"JAX backend: {jax.default_backend()}")

    base_dir = Path(__file__).resolve().parent
    sim_path = base_dir / "simulated_data" / f"sim_data_{dataset_index:03d}.rds"
    params_path = base_dir / "simulated_data" / "true_params.rds"
    out_path = base_dir / "coverage_results" / f"profile_{PROFILE_NAME}_{dataset_index:03d}.rds"

    sim_data = read_sim_data_rds(sim_path)
    true_params = read_params_rds(params_path)
    pomp_dict = build_pomp_dict(sim_data, true_params)
    parameter_shared = generate_parameter_profile(params_path, dataset_index)

    print("Round 1 starting...")
    round_one = run_mif_round(
        pomp_dict=pomp_dict,
        starts=parameter_shared,
        rw_sd_dict=generate_sd(0.05, PROFILE_NAME),
        nmif=200,
        mif_particles=ALGORITHMIC_PARAMS["Mp"][RUN_LEVEL - 1],
        pfilter_particles=ALGORITHMIC_PARAMS["Np"][RUN_LEVEL - 1],
        pfilter_reps=ALGORITHMIC_PARAMS["Np_rep"][RUN_LEVEL - 1],
        key_seed=10_000 + dataset_index,
    )
    print("Round 1 done.")

    round_two_starts = select_round_two_starts(round_one)

    print("Round 2 starting...")
    round_two = run_mif_round(
        pomp_dict=pomp_dict,
        starts=round_two_starts,
        rw_sd_dict=generate_sd(0.04, PROFILE_NAME),
        nmif=300,
        mif_particles=ALGORITHMIC_PARAMS["Mp"][RUN_LEVEL - 1],
        pfilter_particles=ALGORITHMIC_PARAMS["Np"][RUN_LEVEL - 1],
        pfilter_reps=ALGORITHMIC_PARAMS["Np_rep"][RUN_LEVEL - 1],
        key_seed=20_000 + dataset_index,
    )
    print("Round 2 done.")

    output_profile = make_output_profile(round_two)
    write_profile_rds(output_profile, out_path)

    print(f"Saved subset_data_{PROFILE_NAME} for dataset {dataset_index}")
    print(f"  Rows: {len(output_profile)}")
    print(f"  Path: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
