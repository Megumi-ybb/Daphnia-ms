#!/usr/bin/env python3
"""
benchmark_f32.py
================
Test whether float32 is (a) FASTER and (b) FAITHFUL vs float64 for the coverage
filter, before risking it on real outputs.

Precision is chosen by PYPOMP_X64 (1=f64, 0=f32), read at import by the profile
module.  Run it TWICE on the same FREE GPU:

    PYPOMP_X64=1 CUDA_VISIBLE_DEVICES=7 python benchmark_f32.py    # f64 baseline FIRST
    PYPOMP_X64=0 CUDA_VISIBLE_DEVICES=7 python benchmark_f32.py    # then f32 -> prints VERDICT

Each run measures, on simulated dataset b=1:
  * FIDELITY : pfilter log-likelihood across a 20-point `rn` grid (true*0.1 .. true*10,
               other params at truth).  The high-rn end pushes the F/P states toward
               their large magnitudes -- exactly where float32 would lose precision.
  * SPEED    : s/iter of a short mif (M=6) at production width (1800 fits, J=500).

It caches each precision's result next to the script; once both exist it prints a
side-by-side verdict:
   - f32 FAITHFUL if its loglik matches f64 within ~3x the pfilter Monte-Carlo se at
     every grid point (and no NaN/Inf);
   - speedup = f64_s_per_iter / f32_s_per_iter.

Reads only the committed sim CSVs; writes only tiny json caches (_bench_f32/f64.json).
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

import pypomp_coverage_profile as P   # sets jax_enable_x64 from PYPOMP_X64 at import
import jax
import pypomp as pp

PREC = "f64" if jax.config.jax_enable_x64 else "f32"
base = Path(__file__).resolve().parent
sim = P.read_sim_data(base / "simulated_data", 1)
tp = P.read_params(base / "simulated_data")
print(f"=== precision={PREC}  x64={jax.config.jax_enable_x64}  device={jax.devices()[0]} ===",
      flush=True)

# ---------- FIDELITY: pfilter loglik across a 20-point rn grid (others at truth) ----------
NG = 20
grid = np.exp(np.linspace(np.log(tp["rn"] * 0.1), np.log(tp["rn"] * 10.0), NG))
rows = []
for g in grid:
    r = {n: tp[n] for n in P.PARAM_NAMES}
    r["rn"] = float(g)
    rows.append(r)
fid_design = pd.DataFrame(rows)[P.PARAM_NAMES]

panel = pp.PanelPomp(Pomp_dict=P.build_pomp_dict(sim, tp),
                     theta=P.make_panel_parameters(fid_design))
panel.pfilter(J=1000, reps=10, key=jax.random.key(7), chunk_size=1)
ll = np.asarray(panel.results_history[-1].logLiks.values)        # (NG, U, reps)
loglik = np.nansum(P.logmeanexp(ll, axis=2), axis=1)             # (NG,) panel loglik per grid pt
per_rep_panel = np.nansum(ll, axis=1)                            # (NG, reps) sum over units
se = np.nanstd(per_rep_panel, axis=1) / np.sqrt(per_rep_panel.shape[-1])  # (NG,)
n_bad = int(np.sum(~np.isfinite(loglik)))
print(f"fidelity: pfilter loglik over {NG} rn-grid points; "
      f"range [{np.nanmin(loglik):.1f}, {np.nanmax(loglik):.1f}]  non-finite={n_bad}", flush=True)

# ---------- SPEED: short mif s/iter at production width ----------
design = P.generate_parameter_profile(tp, 1, "rn", 60, 30)       # 1800 fits


def mif_once():
    panel_m = pp.PanelPomp(Pomp_dict=P.build_pomp_dict(sim, tp),
                           theta=P.make_panel_parameters(design))
    rw = pp.RWSigma(sigmas=P.generate_sd(0.05, "rn"),
                    init_names=[]).geometric_cooling(a=P.COOLING_A)
    t0 = time.time()
    panel_m.mif(J=500, M=6, rw_sd=rw, key=jax.random.key(42), block=True, vmap_chunk_size=None)
    return time.time() - t0


mif_once()                       # warm (compile)
s_iter = mif_once() / 6.0
print(f"speed: {s_iter:.2f} s/iter", flush=True)

# ---------- cache + verdict ----------
res = {"prec": PREC, "loglik": loglik.tolist(), "se": se.tolist(),
       "s_iter": s_iter, "n_bad": n_bad}
(base / f"_bench_{PREC}.json").write_text(json.dumps(res))
other = base / f"_bench_{'f32' if PREC == 'f64' else 'f64'}.json"
if not other.exists():
    print(f"\n(saved {PREC}; now run the OTHER precision to get the verdict)")
else:
    o = json.loads(other.read_text())
    f64 = res if PREC == "f64" else o
    f32 = res if PREC == "f32" else o
    a = np.array(f64["loglik"]); b = np.array(f32["loglik"]); s = np.array(f64["se"])
    diff = np.abs(b - a)
    tol = 3.0 * np.maximum(s, 1e-6)
    finite = np.isfinite(a) & np.isfinite(b)
    worst = float(np.nanmax(diff[finite])) if finite.any() else float("inf")
    faithful = bool(finite.all() and np.all(diff[finite] <= tol[finite]))
    print("\n==================== VERDICT ====================")
    print(f"  grid points        : {len(a)}   f32 non-finite: {f32['n_bad']}")
    print(f"  max |f32 - f64| ll  : {worst:.3f}   (tol ~3*se median {np.median(tol):.3f})")
    print(f"  FIDELITY           : {'OK  - f32 reproduces the f64 likelihood' if faithful else 'FAIL - f32 drifts the likelihood; keep f64'}")
    print(f"  SPEED              : f64 {f64['s_iter']:.2f}s/iter  f32 {f32['s_iter']:.2f}s/iter"
          f"  -> x{f64['s_iter']/max(f32['s_iter'],1e-9):.2f}")
    if faithful:
        print("  => float32 is SAFE. Re-run the sweep with PYPOMP_X64=0 for the speedup.")
    else:
        print("  => float32 BREAKS fidelity here. Stay on float64 (PYPOMP_X64=1).")
    print("=================================================")
