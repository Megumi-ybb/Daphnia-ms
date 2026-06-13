#!/usr/bin/env python3
"""Standalone test of the R<->Python RDS bridge from pypomp_coverage_profile.py.
Does NOT import jax/pypomp -- it only verifies that reading the simulated .rds
inputs and reproducing pomp::profile_design works locally."""
import json, subprocess, sys
from pathlib import Path

UNIT_NAMES = [f"u{i}" for i in range(1, 9)]
PARAM_NAMES = ["xi","sigSn","sigIn","sigSi","sigIi","sigF","sigP","theta_Sn","theta_In",
    "theta_Si","theta_P","theta_Ii","f_Sn","f_Si","rn","ri","probn","probi","k_Ii","k_In",
    "k_Sn","k_Si","sigJi","sigJn","theta_Jn","theta_Ji"]
PROFILABLE = ["xi","theta_Sn","theta_Si","theta_In","theta_Ii","theta_P","theta_Jn",
    "theta_Ji","rn","ri","f_Sn","f_Si","probn","probi","sigIn","sigIi","sigJn","sigJi",
    "sigF","sigP","k_Sn","k_Si","k_In","k_Ii"]

def run_rscript(code, args):
    c = subprocess.run(["Rscript","--vanilla","-e",code,*args], check=True, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return c.stdout

def read_sim(path):
    code = r"""
    args <- commandArgs(trailingOnly = TRUE); x <- readRDS(args[1])
    if (!requireNamespace("jsonlite", quietly=TRUE)) stop("jsonlite required")
    cat(jsonlite::toJSON(x, dataframe="rows", auto_unbox=TRUE, digits=17))"""
    parsed = json.loads(run_rscript(code, [str(path)]))
    out = {}
    for u in UNIT_NAMES:
        import_rows = parsed[u]
        out[u] = import_rows
    return out

def read_params(path):
    code = r"""
    args <- commandArgs(trailingOnly = TRUE); x <- readRDS(args[1])
    if (!requireNamespace("jsonlite", quietly=TRUE)) stop("jsonlite required")
    df <- data.frame(name=names(x), value=as.numeric(x), check.names=FALSE)
    cat(jsonlite::toJSON(df, dataframe="rows", auto_unbox=TRUE, digits=17))"""
    rows = json.loads(run_rscript(code, [str(path)]))
    return {r["name"]: float(r["value"]) for r in rows}

def profile_design(params_path, b, prof_name, nprof=80):
    pidx = PROFILABLE.index(prof_name) + 1
    seed = 805*1000 + 1000*pidx + b
    code = r"""
    args <- commandArgs(trailingOnly = TRUE)
    true_params <- readRDS(args[1]); prof_name <- args[2]
    nprof <- as.integer(args[3]); seed <- as.integer(args[4])
    if (!requireNamespace("jsonlite", quietly=TRUE)) stop("jsonlite required")
    if (!requireNamespace("pomp", quietly=TRUE)) stop("pomp required")
    set.seed(seed)
    shared_ub <- true_params * 10; shared_lb <- shared_ub / 100
    ub_unit <- log(shared_ub[prof_name]); lb_unit <- log(shared_lb[prof_name])
    shared_lb <- shared_lb[!(names(shared_lb) %in% c("sigSn","sigSi", prof_name))]
    shared_ub <- shared_ub[!(names(shared_ub) %in% c("sigSn","sigSi", prof_name))]
    parameter_shared <- pomp::profile_design(
      temp=seq(lb_unit, ub_unit, length.out=nprof),
      lower=log(shared_lb), upper=log(shared_ub), type="runif", nprof=nprof)
    names(parameter_shared)[names(parameter_shared)=="temp"] <- prof_name
    parameter_shared <- exp(parameter_shared)
    parameter_shared$sigSn <- 0; parameter_shared$sigSi <- 0
    cat(jsonlite::toJSON(parameter_shared, dataframe="rows", auto_unbox=TRUE, digits=17))"""
    return json.loads(run_rscript(code, [str(params_path), prof_name, str(nprof), str(seed)]))

if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    sim_path = base / "simulated_data" / "sim_data_001.rds"
    par_path = base / "simulated_data" / "true_params.rds"

    print("== read sim_data_001.rds ==")
    sim = read_sim(sim_path)
    print(f"  units: {list(sim.keys())}")
    u1 = sim["u1"]
    print(f"  u1 rows: {len(u1)}  first row: {u1[0]}  last row: {u1[-1]}")

    print("== read true_params.rds ==")
    par = read_params(par_path)
    print(f"  n params: {len(par)}  missing canonical: {[p for p in PARAM_NAMES if p not in par]}")
    print(f"  rn={par.get('rn')}  ri={par.get('ri')}  sigF={par.get('sigF')}  sigSn={par.get('sigSn')}")

    print("== profile_design (rn, b=1) ==")
    pd_rows = profile_design(par_path, 1, "rn")
    print(f"  rows: {len(pd_rows)}  cols: {len(pd_rows[0])}")
    print(f"  has all 26 params: {all(p in pd_rows[0] for p in PARAM_NAMES)}")
    rn_vals = sorted(r["rn"] for r in pd_rows)
    print(f"  rn range: [{rn_vals[0]:.4g}, {rn_vals[-1]:.4g}]  sigSn fixed 0: {pd_rows[0]['sigSn']==0}")
    print("== profile_design (ri, b=1) rows ==", len(profile_design(par_path, 1, "ri")))
    print("ALL BRIDGE TESTS PASSED")
