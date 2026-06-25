#!/usr/bin/env python
# ============================================================================
# GPU / pypomp fit for ONE (truth, b, model) of the AIC selection-rate study.
# ----------------------------------------------------------------------------
# Thin wrapper around the VALIDATED bootstrap-LRT GPU port
#   ../bootstrap_lrt/pypomp_bootstrap_lrt.py
# (its model + 2-round mif2 already matches the R fit to ll -880.44 vs -880.56).
# This wrapper changes only THREE things needed by the AIC study, and reuses the
# port's model, mif2 rounds, and result writer verbatim:
#
#   1. SIM DIR by truth.  --truth shared        -> ../coverage_study/simulated_data
#                         --truth unit_specific -> ./simulated_data_unit_specific
#      (both already hold sim_data_all.csv + true_params.csv; see the R converters
#       in pypomp_bootstrap_lrt.make_sim_data_all_csv / the aic_power CSV step.)
#
#   2. ALT TARGET = theta_In ALONE.  The LRT port only ships theta_Ii_theta_In
#      (a 2-param joint block). The AIC study's alt makes theta_In unit-specific
#      on its own -> 8 unit-specific entries, k_alt = 31. We register that target.
#
#   3. AIC.  k_null = 24, k_alt = 31 (project convention: boundary-fixed sigSn,
#      sigSi excluded; +7 for the alt). AIC = 2k - 2*ll. Result saved to
#      results/fit_<model>_<truth>_<b>.rds in collect_aic.R's schema (via pyreadr).
#
# Usage (one fit):
#   PYPOMP_USE_GPU=1 CUDA_VISIBLE_DEVICES=5 \
#     python pypomp_aic.py <b 1..50> <model null|alt> --truth <shared|unit_specific>
#
# The driver run_aic_gpu.sh sweeps all 200 (truth,b,model) with file-skip resume.
# ============================================================================
from __future__ import annotations
import argparse
import os
import sys
from pathlib import Path

# Make the bootstrap_lrt port importable and let it pick the sim dir from env.
_HERE = Path(__file__).resolve().parent
_LRT_DIR = _HERE.parent / "bootstrap_lrt"
sys.path.insert(0, str(_LRT_DIR))


def _sim_dir_for(truth: str) -> Path:
    if truth == "shared":
        return _HERE.parent / "coverage_study" / "simulated_data"
    if truth == "unit_specific":
        return _HERE / "simulated_data_unit_specific"
    raise ValueError(f"truth must be 'shared' or 'unit_specific', got {truth!r}")


def main(argv) -> int:
    ap = argparse.ArgumentParser(description="GPU AIC-study fit for one (truth,b,model).")
    ap.add_argument("b", type=int, help="panel index 1..50")
    ap.add_argument("model", choices=["null", "alt"], help="null (all-shared) or alt (theta_In unit-specific)")
    ap.add_argument("--truth", required=True, choices=["shared", "unit_specific"])
    ap.add_argument("--run-level", type=int, default=int(os.environ.get("PYPOMP_RUN_LEVEL", "3")),
                    choices=[1, 2, 3])
    ap.add_argument("--n-starts", type=int,
                    default=int(os.environ["PYPOMP_NSTARTS"]) if os.environ.get("PYPOMP_NSTARTS") else None)
    args = ap.parse_args(argv)
    if not (1 <= args.b <= 50):
        raise ValueError("b must be 1..50")

    # ---- Resume guard (mirror run_sweep.sh's file-skip) ----
    out_dir = _HERE / "results"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / f"fit_{args.model}_{args.truth}_{args.b}.rds"
    if out_path.exists() and os.environ.get("PYPOMP_OVERWRITE", "0").lower() not in {"1", "true", "yes"}:
        print(f"skip (exists): {out_path}")
        return 0

    # ---- Point the LRT port at the right sim dir BEFORE importing it ----
    sim_dir = _sim_dir_for(args.truth)
    for fn in ("sim_data_all.csv", "true_params.csv"):
        if not (sim_dir / fn).exists():
            raise FileNotFoundError(
                f"{sim_dir/fn} missing. Convert the .rds to CSV first "
                f"(see pypomp_bootstrap_lrt.make_sim_data_all_csv docstring / aic_power CSV step).")
    os.environ["PYPOMP_SIM_DIR"] = str(sim_dir)

    import pypomp_bootstrap_lrt as lrt  # noqa: E402  (after env is set; runs JAX device init)

    # Override the port's hardcoded _sim_dir() with our env-selected directory.
    lrt._sim_dir = lambda: sim_dir  # type: ignore[attr-defined]

    # Register the AIC study's single-param alt: theta_In alone -> 8 unit-specific.
    if "theta_In" not in lrt.ALT_PARAMS_MAP:
        lrt.ALT_PARAMS_MAP["theta_In"] = ["theta_In"]
        lrt.ALT_ORDER = list(lrt.ALT_PARAMS_MAP.keys())

    target = "null" if args.model == "null" else "theta_In"
    rl = args.run_level
    Mp = lrt.ALGO["Mp"][rl]; Np = lrt.ALGO["Np"][rl]; Np_rep = lrt.ALGO["Np_rep"][rl]
    nmif_rounds = lrt.ALGO["Nmif"][rl]
    nstarts = args.n_starts if args.n_starts is not None else lrt.ALGO["nstarts"][rl]
    vmap_chunk = os.environ.get("PYPOMP_VMAP_CHUNK")
    vmap_chunk = int(vmap_chunk) if vmap_chunk else None

    import jax  # noqa: E402
    import pandas as pd  # noqa: E402

    is_null = (target == "null")
    specific_names = [] if is_null else lrt.ALT_PARAMS_MAP[target]
    # Stable seed offset: distinct per truth and per model (shared 0 / unit 500, +100000 alt),
    # matching the R fit scripts' seed structure.
    truth_off = 0 if args.truth == "shared" else 500
    tgt_idx = (0 if is_null else 1) * 1 + (truth_off // 500) * 7  # disjoint key streams

    print(f"=== pypomp AIC fit  truth={args.truth}  b={args.b}  model={args.model} (target={target}) ===")
    print(f"JAX backend = {jax.default_backend()}  devices = {jax.devices()}")
    print(f"sim_dir = {sim_dir}")
    print(f"specific block = {specific_names if specific_names else '(none; all-shared null)'}")
    print(f"Mp={Mp} Np={Np} Np_rep={Np_rep} Nmif={nmif_rounds} nstarts={nstarts}")

    sim_data = lrt.read_sim_data(args.b)
    true_params = lrt.read_true_params()
    theta0 = {n: true_params[n] for n in lrt.PARAM_NAMES}
    pomp_dict = lrt.build_pomp_dict(sim_data, theta0)
    starts = pd.DataFrame([theta0] * nstarts)[lrt.PARAM_NAMES]

    print(f"Round 1 (Nmif={nmif_rounds[0]})...", flush=True)
    round_one = lrt.run_mif_round(
        pomp_dict, starts, specific_names, lrt.make_rw_sd(lrt.RW_ROUNDS[0]),
        nmif=nmif_rounds[0], Mp=Mp, Np=Np, Np_rep=Np_rep,
        key_seed=802_000 + 1000 * tgt_idx + truth_off + args.b, vmap_chunk_size=vmap_chunk)
    starts_df, carry = lrt.select_round_two(round_one, specific_names)

    print(f"Round 2 (Nmif={nmif_rounds[1]})...", flush=True)
    round_two = lrt.run_mif_round_carry(
        pomp_dict, starts_df, carry, specific_names, lrt.make_rw_sd(lrt.RW_ROUNDS[1]),
        nmif=nmif_rounds[1], Mp=Mp, Np=Np, Np_rep=Np_rep,
        key_seed=802_000 + 100_000 + 1000 * tgt_idx + truth_off + args.b, vmap_chunk_size=vmap_chunk)

    best = max(range(len(round_two)), key=lambda i: round_two[i]["loglik"])
    ll = round_two[best]["loglik"]; se = round_two[best]["se"]; coef = round_two[best]["coef"]

    # ---- AIC (project convention: null k=24 excl. boundary sigSn/sigSi; alt +7) ----
    k = 24 if is_null else 31
    aic = 2 * k - 2 * ll
    n_specific = len([c for c in coef if "[" in c])
    if not is_null and n_specific != 8:
        print(f"WARNING: alt has {n_specific} unit-specific entries, expected 8")

    # ---- Save in collect_aic.R's schema via pyreadr ----
    import pyreadr  # noqa: E402
    rec = {"b": [args.b], "truth": [args.truth], "model": [args.model],
           "ll": [ll], "se": [se], "k": [k], "AIC": [aic]}
    if not is_null:
        rec["alt_name"] = ["theta_In"]
    for cname, cval in coef.items():
        rec[f"coef.{cname}"] = [float(cval)]
    pyreadr.write_rds(str(out_path), pd.DataFrame(rec))

    print("\n================= RESULT =================")
    print(f"truth={args.truth} b={args.b} model={args.model}")
    print(f"ll = {ll:.4f}  se = {se:.4f}  k = {k}  AIC = {aic:.4f}")
    print(f"coef entries = {len(coef)} (unit-specific {n_specific})")
    print(f"saved -> {out_path}")
    print("==========================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
