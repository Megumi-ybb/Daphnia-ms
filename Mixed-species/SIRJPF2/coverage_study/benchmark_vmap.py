#!/usr/bin/env python3
"""
benchmark_vmap.py
=================
Measure whether PYPOMP_VMAP_CHUNK (running the 8 panel units in PARALLEL via vmap
instead of the default SEQUENTIAL lax.scan over units) makes each mif iteration
cheaper on the GPU -- without committing the whole 4-day sweep.

It builds a panel at PRODUCTION width (default 1800 fits, J=500, 8 units) and times
a short mif (default M=6 iterations) twice per setting (1st warms the compile, 2nd is
timed), for vmap_chunk_size in {None (default, sequential), 8 (parallel)}.

Run on a FREE GPU (takes ~15-25 min):
    CUDA_VISIBLE_DEVICES=7 python benchmark_vmap.py
Faster/rougher:
    BENCH_NSTARTS=15 BENCH_M=4 CUDA_VISIBLE_DEVICES=7 python benchmark_vmap.py

Reads only the committed sim CSVs; writes nothing. Importing the profile module
configures JAX (float64 + GPU) exactly as the real run does.
"""
import os
import time
from pathlib import Path

# importing the profile module sets jax_enable_x64 + GPU BEFORE jax import, and
# gives us the exact model-construction helpers the real sweep uses.
import pypomp_coverage_profile as P  # noqa: E402
import jax  # noqa: E402
import pypomp as pp  # noqa: E402

NPROF = 60
NSTARTS = int(os.environ.get("BENCH_NSTARTS", "30"))   # 60*30 = 1800 = production width
M = int(os.environ.get("BENCH_M", "6"))                # short mif just for timing
J = P.ALGORITHMIC_PARAMS["Mp"][2]                      # 500 (run-level 3)
PARAM, B = "rn", 1

base = Path(__file__).resolve().parent
sim = P.read_sim_data(base / "simulated_data", B)
tp = P.read_params(base / "simulated_data")
design = P.generate_parameter_profile(tp, B, PARAM, NPROF, NSTARTS)

print(f"device   : {jax.devices()[0]}   x64={jax.config.jax_enable_x64}")
print(f"benchmark: {len(design)} fits, J={J}, M={M}, 8 units\n")


def run_once(vchunk):
    pomp_dict = P.build_pomp_dict(sim, tp)             # fresh panel (mif fits in place)
    panel = pp.PanelPomp(Pomp_dict=pomp_dict, theta=P.make_panel_parameters(design))
    rw = pp.RWSigma(sigmas=P.generate_sd(0.05, PARAM),
                    init_names=[]).geometric_cooling(a=P.COOLING_A)
    t0 = time.time()
    panel.mif(J=J, M=M, rw_sd=rw, key=jax.random.key(42),
              block=True, vmap_chunk_size=vchunk)
    return time.time() - t0


def timed(vchunk, label):
    print(f"[{label}] warming (compile)...", flush=True)
    try:
        run_once(vchunk)                               # warm: pays the compile
        t = run_once(vchunk)                           # timed: cache hit
    except Exception as e:
        print(f"[{label}] FAILED: {type(e).__name__}: {e}\n")
        return None
    print(f"[{label}] {t:.1f}s total  ->  {t / M:.2f} s/iter\n", flush=True)
    return t / M


a = timed(None, "default  units SEQUENTIAL (vmap_chunk=None)")
b = timed(8, "vmap_chunk=8  units PARALLEL")

print("=" * 56)
if a and b:
    print(f"per-iter: default {a:.2f}s   chunk8 {b:.2f}s   ->  x{a / b:.2f} "
          f"({'FASTER' if b < a else 'no gain'})")
    print("If FASTER: validate one FULL job both ways (compare output max loglik)")
    print("before rolling PYPOMP_VMAP_CHUNK=8 to the sweep.")
else:
    print("one setting failed - see above (chunk=8 may not support this config)")
print("=" * 56)
