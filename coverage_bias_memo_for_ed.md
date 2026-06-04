Subject: SIRJPF2 coverage — bias diagnostic confirms it's a finite-sample property, not particles

Ed,

Following up on your question of whether the MCAP undercoverage reflects a true property of the model / finite sample or a finite-particle artifact. I ran the three-log-likelihood bias diagnostic on five representative bootstrap datasets (b in {2, 25, 50, 75, 95}). The answer is firm where it counts: the undercoverage is finite-sample displacement of the MLE, not an artifact of under-converged profiling — though I'm careful below not to overclaim that more particles "cannot" help (that would need a particle-count ablation I haven't run).

Two lines of evidence:

1. Unconstrained search beats the truth. I warm-started a 20-chain mif2 run at theta_true with all 24 parameters free. The best fit beats ll(truth) by mean +4.47 log-likelihood (per-dataset best gains: ds2 +5.11, ds25 +3.54, ds50 +6.46, ds75 -0.13, ds95 +7.38). The selection-robust version of this: on 2 of the 5 datasets (ds50, ds95) *every one* of the 20 chains independently exceeds the truth (median chain gain +2.0 and +4.9), so the MLE is genuinely above truth there, before MCAP enters. On ds2/ds25 the "best" gain is largely a max-over-chains selection effect, and ds75 stayed at truth. So treat +4.47 as an upper-leaning (max-over-chains) estimate of E[ll_MLE - ll_truth], not an unbiased one.

2. The coverage profiles are at least as thorough as that global search. Across all 5 datasets x 4 target params, profile_max - ll(truth) runs +5.2 to +12.5 (mean ~+8.5), never negative; and ll_global_best - ll_profile_max is about -4 (negative): the production profiles out-search the 20-chain warm-from-truth diagnostic. So the undercoverage is NOT an artifact of under-converged profiling relative to that baseline. (I am careful here: this does not prove either search reached the true global MLE — profile_max > global_best is partly just a broader max — so I would not claim "more particles cannot help" without a particle/iteration ablation.)

Conclusion: the undercoverage is finite-sample displacement of the MLE from the truth on weakly-identified ridges. The two gaps (~4.5 from the 20-chain search, ~8.5 from the broader profile) are searches of differing thoroughness that bracket the same quantity, ll_MLE - ll_truth. Both sit below the heuristic asymptotic envelope p/2 = 12 (p = 24 free params) — consistent with, but not a derived bound for, finite-sample bias in a weakly-identified model. (Both are maxima over Monte-Carlo replication noise (~0.5 units) and, more substantially, over chains/grid points, so they upward-bias the gap; MCAP corrects the MC component in the interval.) The practical upshot stands: the displacement is real, it is concentrated in weakly-identified directions, and it is not a numerical artifact of under-converged profiling.

Where the bias lives: well-identified params drift modestly (mean log-scale, SD): rn +1.11 (0.92), probn +0.92 (0.76), f_Sn -0.93 (0.77). Weakly-identified / ridge params drift hard: ri -1.18 (2.56), probi -1.25 (2.89), f_Si +1.41 (2.73), sigIn +2.08 (4.51), sigIi +3.66 (6.89). Within-dataset chain SD tells the same story: rn ~1.5 vs sigIn/sigIi ~4.8.

One wrinkle worth your eye: the picture is not perfectly uniform. r_n fits the displacement story cleanly (MLE biased high, CIs miss). sigma_F, by contrast, has near-zero diagnostic drift yet still undercovers (74%), which points to interval width rather than MLE displacement for that parameter — useful to know, since a width problem is exactly the kind of thing a cutoff calibration (option B) could correct if we decide we want one.

My recommendation: lead with C — report the honest empirical coverage (with Monte Carlo SE) under the exploratory-not-confirmatory framing you suggested — and keep cutoff calibration (B) only as a sensitivity check, not as the headline. Reasons:

- It matches what the referee actually asked for ("an exploration of these issues") and your own read that MCAP is doing well apart from a small bias.
- B=100 is on the small side to calibrate a 95% quantile reliably, so honest coverage + SE is the firmer claim; calibration is better shown as "the undercoverage is correctable if desired" than as the reported intervals.
- Calibration can't rescue r^i (non-identifiable product ridge), so a blanket "we calibrated and everything now covers" wouldn't be truthful; C handles r^i honestly.
- It avoids rewriting the Table 1 main-text CIs — I'd leave those as nominal MCAP, with the S3 pointer carrying the caveat.

So unless you'd rather make calibration the headline, I'll write S3 that way: honest coverage + exploratory framing as primary, bootstrap-calibrated cutoff as a sensitivity check, Table 1 unchanged. The one judgment call I'd value your read on: whether the sensitivity-check calibration is even worth showing given the B=100 noise, or whether we report honest coverage and stop there.

Bo
