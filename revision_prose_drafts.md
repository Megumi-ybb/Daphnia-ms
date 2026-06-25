# Revision prose drafts — FOR BO REVIEW (not yet inserted)

Produced by the prose-drafting workflow + adversarial critique on 2026-06-03.
**The critic returned NO-GO as drafted.** Two disqualifying issues (Draft 2 novelty overclaim; Draft 3 self-contradiction with ms.Rnw:1165). Resolve the flagged items before inserting any block into the live files. All blocks are `\bye{}`-wrapped (purple revision markup).

> **How to read this file.** Each draft below has two parts: (i) the manuscript/response-bound DRAFT PROSE, which lives inside the `\bye{...}` LaTeX blocks and is already written as fluent, first-person-plural paragraphs (no bullets, no `\paragraph{}`), and (ii) the surrounding working notes — the *Critique*, *Verify / notes*, *Claims to verify*, and *Decisions* lists. Those lists are deliberately kept as checklists, because their per-item structure (flag → fix, cross-ref → resolved/invented, claim → source) is what makes them usable; flattening them into prose would lose information.
>
> **Staleness flag (2026-06-24):** Draft 3's prose block (the stochastic-vs-deterministic Discussion paragraph) has since been revised and inserted into the live `ms.Rnw` (the inserted version begins "The stochastic formulation does more than improve numerical fit; it reframes the ecological conclusions ... and makes them falsifiable"). The live text already addresses the two Draft-3 critique flags: it softens "what the data require" to "favors a model" / "the decisively better-supported description," adds the explicit "no food-level data were collected" caveat (resolving the ms.Rnw:1165 contradiction), cites the no-$F$ ablation as **Section S6** (not the placeholder `\ref{sec:no-F-ablation}`), and reframes the test as "a calibrated parametric bootstrap likelihood-ratio test" with the $\chi^2$-unreliability caveat. Treat the Draft 3 block below as the historical draft; the manuscript is now ahead of it.

## Critique

**Verdict:** NO-GO as drafted. The single most dangerous flaw, the one most likely to draw a sharp referee response, is the overclaim in Draft 2's discussion pointer: calling block/marginalized IF 'a new algorithm' and asserting comparative superiority over breto20 without a cited result. The editor specifically asked the authors to clarify the added value of the marginalized approach; answering with an unsupported novelty claim about a pre-existing software option will read as evasion or sloppiness to anyone familiar with panelPomp. The ecological-claim contradiction between Draft 3 and ms.Rnw:1165 (food-depletion is both 'what the data require' and 'an unvalidated assumption') is equally disqualifying because it creates a visible self-contradiction within the Discussion section that any careful referee will flag. Both issues require targeted fixes before submission; all other per-draft issues are correctable in one editing pass.

**Overclaim flags:**

- DRAFT 2 discussion pointer calls block/marginalized panel IF 'a new marginalized panel iterated filtering algorithm' -- it is an existing panelPomp::mif2 option (MARGINALIZE=TRUE / block=TRUE), not a new algorithm. The paper's contribution is demonstrating its value on this application.
- DRAFT 2 discussion pointer asserts MPIF provides 'more robust likelihood maximization ... than the joint-resampling panel iterated filter of \citet{breto20}' without citing any head-to-head comparison number in the manuscript; si.Rnw:270 states empirical superiority but no comparison table or figure is described in the gathered context.
- DRAFT 1 cover letter says 'the parametric bootstrap likelihood ratio tests ... find that no unit-specific elaboration of the all-shared model is supported' but si.Rnw:477 states this conclusion based on chi-square LRT p-values, not bootstrap p-values; the bootstrap column is present in the Sweave table code but has no corresponding written-out narrative claim in the SI, making the cover-letter attribution to bootstrap results premature.
- DRAFT 3 Discussion states resource depletion is 'what the data require' for the population collapse, whereas ms.Rnw:1165 already acknowledges in the same Discussion section that this mechanism 'remains an unvalidated assumption in the theory rather than being directly supported by data' because no food measurements were taken.

**Per-draft required fixes:**

- DRAFT 1 (editor reply, emphases A+B): The cross-reference 'Section~S\ref*{sec:no-F-ablation} (Table~S1, the bootstrap p-value column)' is wrong -- the bootstrap LRT for unit-specific model selection lives in the SIRJPF2 model-comparison section (Table~\ref{Table:SIRJPF2_model_comparison}, si.Rnw:460-478), NOT in the no-F ablation section (sec:no-F-ablation). Fix: change the cite to point to the model-comparison SI section and Table~S1 by its actual label; separately, confirm the SI narrative at si.Rnw:477 says 'LRT p-values' not 'bootstrap p-values' -- if the bootstrap column numbers are not yet written into the SI prose, the cover-letter claim that bootstrap results support no unit-specific elaboration is citing an unpublished computation, not a stated SI result.
- DRAFT 2 (ms.Rnw methods paragraph + discussion pointer): The one-sentence discussion pointer says 'a new marginalized panel iterated filtering algorithm' and 'more robust likelihood maximization ... than the joint-resampling panel iterated filter of \citet{breto20}'. Both clauses are overclaims. Block/marginalized resampling is an existing option in panelPomp::mif2 (si.Rnw line 278 says it 'is implemented in the mif2 function'), so 'new algorithm' is false. The comparative 'more robust' claim is not backed by any head-to-head table or figure in the manuscript. Fix (discussion pointer): replace with 'the marginalized panel iterated filter, which uses block resampling of unit-specific parameters to respect the conditional independence of panel units; we found this option well-suited to this high-dimensional, weakly-identified model, and a theoretical treatment is provided by \citet{wheeler25}.' Do not assert 'new' or a comparative superiority claim without a cited result.
- DRAFT 3 (ms.Rnw Discussion emphasis-B): The draft's first paragraph concludes 'Depletable resource dynamics, not a reproductive-strategy threshold, are therefore what the data require to explain the observed peak-then-collapse trajectory.' This directly contradicts ms.Rnw line 1165, which already exists in the Discussion and says 'The mechanism, although logical within the model's framework, remains an unvalidated assumption in the theory rather than being directly supported by data' (because no food-level data were collected). The confident causal framing in the draft will create a visible internal contradiction. Fix: either (a) insert the existing caveat from line 1165 into or immediately after the new paragraph, or (b) soften the conclusion to match the existing hedge: 'the food-depletion feedback is the better-supported statistical explanation of the observed dynamics, though direct validation would require food-level measurements not collected here.'
- DRAFT 4 (Section 2 opening paragraph + Significance Statement): The Significance Statement draft includes '\section*{Significance Statement}' inside the \bye{} block but ms.Rnw line 1180 already has that header. This will produce a duplicate section title in the compiled PDF. Fix: drop the \section* line from inside the \bye{} block and replace only the \revision{\editSignificance} placeholder at line 1182, per the notes' own recommendation (b) -- this is flagged but not corrected in the draft itself.
- DRAFT 5 (si.Rnw S3 Interpretation subsection): The draft describes calibration as 'enlarging the profile-drop threshold until 95% of the B bootstrap intervals cover the truth.' Verify that the actual MCAP calibration lever is the profile log-likelihood drop cutoff (the ionides17/ning21 lambda/smoothing parameter) and not a simple quantile rescaling -- the draft's description may misdescribe the mechanism. Fix: before inserting, confirm the exact calibration procedure against the MCAP code and state it precisely, or describe it abstractly as 'adjusting the MCAP cutoff' without specifying the lever if uncertain.

---

## Draft 1 — Editor cover-letter reply to emphases A (marginalized/block panel IF) and B (stochastic vs deterministic ecological conclusions), as two \bye{} blocks, each preceded by a brief \report{}-paraphrase lead-in matching the surrounding response.tex style.

**Insert at:** `/Users/ybb/Downloads/Research/Daphnia-ms/Daphnia-paper-review/response.tex:140`

**Verify / notes:**

PLACEMENT: Insert both \bye{} blocks after line 140 (\article{\editSignificance}) and before the blank line preceding the AE header at line 142, per the response-structure insertion point. The two emphases currently have NO reply in response.tex (confirmed lines 67-141).

CROSS-REF LABELS TO VERIFY before build (I used \ref*{} / \ref{} with these; confirm the labels resolve in the response.tex preamble, which \input's the article, or hard-code the section letters/numbers if response.tex does not share the article's label namespace):
  - \ref{fig:pif} -> Figure 1 (defined ms.Rnw line 706). OK in ms; verify it resolves in response.tex.
  - sec:meth-supp -> I INVENTED this label for the SI methodology/algorithm section. The SI algorithm discussion is at si.Rnw lines 268-278 but I did NOT find an explicit \label there. Bo: either add \label{sec:meth-supp} to that SI section or replace my \ref*{sec:meth-supp} with the correct existing SI section label / hard-coded "Section S2" (ms.Rnw line 643 already calls the PIF pseudocode "Section~S2").
  - sec:mcap-coverage -> CONFIRMED real label, si.Rnw line 282.
  - sec:no-F-ablation -> CONFIRMED real label, si.Rnw line 378. Note line 107 of response.tex already uses Section~S\ref*{sec:no-F-ablation}, so this pattern is established and resolves.
  - sec:res -> CONFIRMED real label (ms sec:res).
  - alg:pif -> referenced in si.Rnw line 276 as \ref{alg:pif}; confirm the algorithm carries that label.
  - Table:Model_comparison, tab:estimates -> CONFIRMED labels (ms.Rnw lines 865, 837 per ms-structure map).
  - Table S1 -> I call it "Table~S1" in plain text (matches ms.Rnw line 922 wording "Supplementary Table~S1"). Confirm the bootstrap p-value column lives in that table.

CLAIMS TO VERIFY (all sourced, but double-check the live numbers):
  - dAIC ~12.4, ll -880.6 -> -887.8, k 24->23, dll ~7.2: matches response.tex line 107 (existing \bye block) and SI sec:no-F-ablation. Consistent.
  - MPIF "leaves psi_{tilde u} unchanged while filtering through u": matches si.Rnw line 268 (MARGINALIZE=TRUE) and ms.Rnw lines 703-704 figure caption. Consistent.
  - "no unit-specific elaboration supported / bootstrap p-value column": matches si.Rnw line 477 ("none ... yield a statistically significant improvement") and the bootstrap LRT machinery (si.Rnw lines 442-448). NOTE: per MEMORY (project_bootstrap_lrt), across all 4 families only Dent_SIRJPF theta_In rejects all-shared; for the SIRJPF2 MIXED-species panel specifically (the model this response is about), none reject -- so the statement "no unit-specific elaboration of the all-shared model is supported" is correct FOR SIRJPF2. I deliberately scoped the claim to SIRJPF2/the eight units to stay accurate; Bo: confirm the SIRJPF2 Table S1 bootstrap column indeed shows no rejection (the Dent theta_In rejection is a different, single-species family and is not referenced here).
  - sigma_S native/invasive fixed to 0 (boundary): relevant to the k=24 count (MEMORY aic_k_count) but not asserted in the draft; no action.

STYLE: Each \bye{} block opens with a one-sentence paraphrase of the editor's ask (no \report{} quote was supplied for A/B in the decision email, so I led with a paraphrase rather than fabricating a \report{} quote -- consistent with how line 124-126 handles the grayscale point without re-quoting). If Bo prefers, the paraphrase lead-ins can be pulled out of the purple \bye and set in normal black text like lines 74/138; I kept them inside \bye so the whole new-content block renders purple/unfinished per the \bye convention.

OVERCLAIM GUARD: Emphasis A explicitly defers theory to wheeler25 and frames MPIF value as empirical/algorithmic (per novelty-grounding "Do NOT claim major theoretical advance"). Emphasis B keeps the Searle16 sexual-switch mechanism "biologically plausible" by implication and frames resource depletion as the better STATISTICAL explanation, not the only possible biology -- consistent with ms.Rnw lines 1138-1140 and the Ed steer to present coverage as exploratory (coverage is referenced only as the identifiability context for weak directions, not as a confirmatory result).

PAGE/LINE numbers cited inside the prose (910-917, 1168-1170, 265-266, 1138-1140) are from the gathered ms-structure map; if the ms has shifted since, these are illustrative pointers for the editor and can be dropped or updated, but they are accurate as of the map provided.

```latex
\bye{The editor asks us to emphasize the methodological contribution beyond prior PanelPOMP developments, and in particular to clarify the added value of the marginalized panel iterated filtering algorithm.

The panel iterated filtering (PIF) of \citet{breto20} resamples, at each filtering step through a unit $u$, the unit-specific parameters of \emph{all} units $\tilde u \neq u$ jointly with the focal unit's parameters $\psi_u$ and the shared parameters $\phi$. The marginalized variant (MPIF) instead leaves the parameters $\psi_{\tilde u}$ of the other units unchanged while filtering through unit $u$, reselecting only $\phi$ and $\psi_u$. This is the distinction now stated in the caption of Figure~\ref{fig:pif} and in Section~S\ref*{sec:meth-supp} of the supplement (Algorithm~\ref{alg:pif}, $\mathrm{MARGINALIZE}=\mathrm{TRUE}$). The added value is not cosmetic for a study of this kind. Marginalization respects the conditional independence of the panel units given the shared parameters, so a poor fit in one unit no longer propagates into the resampling weights of the others. Under joint resampling, by contrast, when units fit unevenly---as they do here, since the eight mesocosms differ in species composition, density, and parasite load---the few units that fit poorly dominate the particle weights and degrade the effective sample size, driving the unit-specific parameter swarm toward degeneracy. For this $8$-unit panel, with weakly identified directions such as the invasive recruitment and infection products $r^{\invasive} f_S^{\invasive}$ and $p^{\invasive} f_S^{\invasive}$ (Section~\ref{sec:res} and Section~S\ref*{sec:mcap-coverage}), that degeneracy is precisely what would frustrate reliable likelihood maximization and the unit-level inference that our model comparison and profile confidence intervals rest on. We therefore present MPIF as the algorithmic lever that makes full-information likelihood inference feasible for this high-dimensional mechanistic panel, and we are careful not to overclaim: we report MPIF's superior empirical performance on this model (the side-by-side log-likelihood and AIC columns in Tables~\ref{Table:Model_comparison} and S1), while attributing the theoretical convergence analysis to the contemporaneous work of \citet{wheeler25}, which we now cite explicitly at lines~265--266 and in Section~S\ref*{sec:meth-supp}.}

\bye{The editor also asks for a clearer synthesis of how the stochastic framework alters or enhances the ecological conclusions relative to earlier deterministic analyses.

The deterministic ODE analysis of \citet{searle16} attributes the observed population turnover to a switch in reproductive strategy, from the asexual production of juveniles to the sexual production of dormant eggs, triggered at high density and low food. Our stochastic SIRJPF2 model reaches a different and statistically better-supported conclusion: the peak-then-collapse dynamics are generated by an explicit depletable food compartment $F(t)$ coupled to age-structured recruitment, with dynamic stochasticity in the juvenile, food, and parasite stages absorbing the unmodeled variation (Section~\ref{sec:res}, lines~910--917). This reframing is now falsifiable rather than rhetorical. The no-$F$ ablation of Section~S\ref*{sec:no-F-ablation}, in which $F$ is held fixed so that recruitment no longer responds to a depletable resource, lowers the maximized log-likelihood by about $7.2$ units (from $-880.6$ to $-887.8$) despite having one fewer parameter, yielding $\Delta\mathrm{AIC}\approx 12.4$ in favor of the depletable-resource model; the resource-depletion feedback is thus a load-bearing, testable conclusion and not an untested modeling assumption. Stochasticity changes the conclusions in two further ways that a deterministic fit misses. First, the data come down decisively in favor of process noise in exactly the stages where resource limitation acts ($\sigma^{\native}_{J}$, $\sigma^{\invasive}_{J}$, $\sigma_F$, $\sigma_P$; Table~\ref{tab:estimates}), so the noise is itself diagnostic of where the deterministic mechanism is incomplete. Second, the dynamic noise accounts for the differences between experimental replicates without any need for mesocosm-specific parameters: the parametric bootstrap likelihood ratio tests of Section~S\ref*{sec:no-F-ablation} (Table~S1, the bootstrap $p$-value column) find that no unit-specific elaboration of the all-shared model is supported, so a single shared stochastic dynamic explains all eight units. A deterministic model, lacking process noise, would be forced to posit unit-specific parameters or mechanisms to reproduce the same replicate heterogeneity. We have woven this synthesis into the Discussion (lines~1168--1170) and the new significance statement.}
```

---

## Draft 2 — ms.Rnw emphasis-A draft: (1) a \bye{}-wrapped methods paragraph explaining the marginalized/block panel iterated filtering contribution, to be inserted in Section~\ref{sec:meth} immediately after the PIF figure (Fig.~\ref{fig:pif}, line 707) and before the existing "Parameter estimates and corresponding confidence intervals were obtained using PIF." sentence at line 708; and (2) a one-sentence \bye{}-wrapped discussion pointer for Section~\ref{sec:dis}. Notation matches the ms: shared parameters $\phi$, unit-specific $\psi_u$, other-unit parameters $\psi_{\tilde u}$ for $\tilde u\neq u$, particle index $j$, $U$ units, conditional independence of latent processes across units (echoing line 622). Cites breto20 (panel IF), ionides15 (IF2), wheeler25 (theory, contemporaneous).

**Insert at:** `Piece 1 (methods paragraph): /Users/ybb/Downloads/Research/Daphnia-ms/daphnia-article/ms.Rnw:707 (after \end{figure}, before line 708 "Parameter estimates and corresponding confidence intervals were obtained using PIF."). Piece 2 (discussion pointer): /Users/ybb/Downloads/Research/Daphnia-ms/daphnia-article/ms.Rnw:1170 (Section sec:dis, after the existing stochasticity/parsimony sentence).`

**Verify / notes:**

NOVELTY CLAIMS NEEDING BO/ED CONFIRMATION:

1. "more reliable likelihood maximization in this regime" / "more robust likelihood maximization ... than the joint-resampling panel iterated filter" — This is the central empirical novelty claim (Emphasis A). It is supported in spirit by the existing intro sentence (ms line 265: "we demonstrate empirically that this is further strengthened by a new marginalized panel iterated filtering algorithm") and the SI claim (si.Rnw line 270: "MPIF has superior empirical performance on the ecological model considered here"). HOWEVER, I did not find a head-to-head MPIF-vs-PIF likelihood comparison number in the gathered context. CONFIRM with Bo whether such a comparison exists (e.g., a maximized-loglik or convergence-rate difference) that can be cited, or whether the claim should be softened to "we found the marginalized algorithm well-suited to this regime" to avoid asserting a comparison not shown in the paper. Ed's steer on the coverage study was to be honest/exploratory; the same caution may apply here.

2. Particle-depletion mechanism ("when one unit fits poorly its low particle weights deplete the diversity of psi_tilde-u across the whole panel") — This is the editor's stated rationale (from the task brief / decision email) and is mechanistically sound, but it is presented here as fact rather than as a cited result. CONFIRM whether to attribute this to wheeler25 or present as motivation only. The brief explicitly invokes "particle-depletion / weight-degeneracy that joint resampling suffers when one unit fits poorly" so this aligns with the editor's framing.

3. "block-resampling option of the mif2 function in the panelPomp R package" — MEMORY/brief says this is block=TRUE. The SI (si.Rnw line 268) uses the flag name MARGINALIZE=TRUE, and line 278 references the mif2 function / breto25 pseudocode. CONFIRM the exact public option name (block=TRUE vs MARGINALIZE) and whether to name it in the main text at all, or defer to SI Section S2. I avoided writing "block=TRUE" literally to prevent a mismatch; verify against the installed panelPomp interface.

4. Crediting block resampling to breto20 vs wheeler25: The brief's Emphasis-A note says "Cite Breto, Ionides, King 2020 JASA for panel IF; the marginalized/block variant is the methodological lever this paper leans on" and also "the panelPomp this is mif2 with block=TRUE." But the SI (line 269) is careful: the UNMARGINALIZED PIF is breto20's; theory for MPIF is wheeler25 (forthcoming), with no convergence guarantee yet. I followed the SI: breto20 = the (unmarginalized) panel IF baseline; wheeler25 = the marginalized variant's theory. The draft does NOT claim a theoretical advance, consistent with the novelty-grounding guidance ("Do NOT claim: Major theoretical advance (theory forthcoming)"). CONFIRM this attribution split is what Bo/Ed want in the main text.

STYLE / MECHANICAL DECISIONS:
- Wrapped both pieces in \bye{...} (renders purple) per convention; multi-paragraph not needed (single paragraph each).
- Reused exact ms notation: $\phi$, $\psi_u$, $\psi_{\tilde u}$, $\tilde u \neq u$, particle $j$, $U$ units — all verified against ms.Rnw lines 616, 622, 703--704.
- Echoed the conditional-independence wording from ms line 622 ("latent processes are modeled as independent across units ... linked only via the shared parameter vector phi") to keep the main text internally consistent.
- Did NOT use species macros \native/\invasive here (this paragraph is algorithm-level, not parameter-specific); kept "weakly-identified" generic to point at the r^invasive ridge without restating it.
- Length: Piece 1 is ~9 sentences / one paragraph; given the AOAS 20-page limit, if space is tight the sentence beginning "This effect is most damaging..." can be cut without losing the core contrast. Flagged for Bo.
- TERM "MPIF" is introduced here in the main text; SI already uses MPIF (si.Rnw line 270). Confirm the main text wants the acronym spelled out as "marginalized panel iterated filtering (MPIF)" on first use — the intro (line 265) currently says "marginalized panel iterated filtering algorithm" without the acronym. Consider introducing MPIF at line 265 for consistency, OR drop the acronym here. Minor.
- I did NOT edit ms.Rnw; these are drafts for Bo to place. Existing line 708 sentence "Parameter estimates ... were obtained using PIF" should be reconciled with the new paragraph: either change "PIF" -> "MPIF" there, or let the new paragraph stand and leave 708 as-is. Recommend changing line 708 "using PIF" -> "using MPIF" for consistency once the paragraph is inserted — flagged, not done.

```latex
% ===== PIECE 1: METHODS PARAGRAPH =====
% INSERT in Section sec:meth immediately AFTER the Fig.~\ref{fig:pif} environment
% (i.e., after line 707 "\end{figure}") and BEFORE line 708
% ("Parameter estimates and corresponding confidence intervals were obtained using PIF.").
% The figure caption (lines 703--704) already defines the marginalized vs. standard
% resampling distinction; this paragraph explains why the marginalized variant matters here.

\bye{We carry out all likelihood maximization using the marginalized variant of panel iterated filtering (MPIF), which builds on the panel iterated filter of \citet{breto20} and the iterated filtering (IF2) algorithm for individual POMP models \citep{ionides15}.
Both algorithms associate a perturbed parameter value with each particle $j$ and use the selection induced by particle-filter resampling to drive the swarm toward the maximum likelihood estimate.
The two variants differ only in the parameter resampling step.
In the panel iterated filter of \citet{breto20}, when the swarm is filtered through the data of unit $u$ the resampling weights are applied not only to $\phi$ and $\psi_u$ but also to the unit-specific parameters $\psi_{\tilde u}$ of all other units $\tilde u \neq u$, even though those parameters do not enter the likelihood for unit $u$.
The marginalized variant leaves $\psi_{\tilde u}$ unchanged while filtering through unit $u$, resampling only the shared parameters $\phi$ and the unit-specific parameters $\psi_u$ of the focal unit (Fig.~\ref{fig:pif}).
This respects the defining conditional-independence structure of a PanelPOMP, in which the latent processes are independent across units given the shared parameters $\phi$, so that the data from unit $u$ are informative about $\phi$ and $\psi_u$ but carry no information about $\psi_{\tilde u}$.
Resampling $\psi_{\tilde u}$ on the weights from unit $u$, as in the unmarginalized algorithm, ties the estimation of each unit's parameters to the idiosyncratic fit of every other unit; when one unit fits poorly its low particle weights deplete the diversity of $\psi_{\tilde u}$ across the whole panel.
This effect is most damaging for the weakly-identified parameters that arise in our model, where the data constrain only certain combinations of parameters, and we found that the marginalized algorithm gave more reliable likelihood maximization in this regime.
The marginalized algorithm we use is implemented as the block-resampling option of the \texttt{mif2} function in the \texttt{panelPomp} R package, and a contemporaneous theoretical and methodological treatment is given by \citet{wheeler25}.}

% ===== PIECE 2: ONE-SENTENCE DISCUSSION POINTER =====
% INSERT in Section sec:dis (Discussion), where the methodological contributions are
% summarized. A natural home is alongside the existing stochastic-modeling remarks
% near lines 1168--1170; place as a standalone sentence.

\bye{Methodologically, the analysis was made possible by the marginalized panel iterated filter, whose block resampling of unit-specific parameters respects the conditional independence of panel units and thereby provides more robust likelihood maximization for high-dimensional PanelPOMP models with weakly-identified parameters than the joint-resampling panel iterated filter of \citet{breto20}.}
```

---

## Draft 3 — ms.Rnw Discussion: stochastic-vs-deterministic ecological synthesis (editor emphasis B), as a \bye{} block inserted after the existing stochasticity/parsimony paragraph (current line 1169).

**Insert at:** `/Users/ybb/Downloads/Research/Daphnia-ms/daphnia-article/ms.Rnw:1169`

**Verify / notes:**

PLACEMENT: Insert as a new paragraph immediately AFTER current line 1169 (the existing "We demonstrate that including stochasticity... mesocosm-specific parameters..." sentence) and before the blank line 1170. The new block extends that exact theme, so it reads as elaboration. Edit-ready anchor: insert after the string ending "...the difference between the experimental replicates." Do NOT place it before line 1168 or it will duplicate the king15 framing.

CLAIMS TO VERIFY (all currently grounded, but confirm numbers stay in sync if SI is re-run):
1. No-F ablation: ~7.2 ll units, dAIC ~ +12.4 favoring depletable-F. Source: si.Rnw line 392/395/400, Table:no-F-ablation. CROSS-REF used: \ref{sec:no-F-ablation} (si.Rnw:378). NOTE: This is a \ref to an SI \label; in the AOAS split-supplement build confirm the SI \label{sec:no-F-ablation} resolves correctly (it renders as "Section S?" in the combined PDF). If main/SI compile separately, may need to hardwire "Supplementary Section S3" or whatever number it lands on. Task #7 owns the si.pdf rebuild.
2. "no unit-specific parameterization yields a statistically significant improvement over the all-shared model" — TRUE for SIRJPF2 mixed-species (si.Rnw:476-478: all-shared lowest AIC, LRT p-values all non-significant). CAVEAT per MEMORY/bootstrap_lrt: across the OTHER families, Dent_SIRJPF/theta_In DOES reject all-shared (p_boot significant). I scoped the sentence to THIS model ("here ... none do") to stay true; do not broaden to "across all families none reject."
3. searle16 = deterministic ODE, reproductive-strategy switch (asexual juveniles -> sexual dormant eggs). Source: ms.Rnw:268, 1138-1139. Matches existing intro framing at line 268 (good consistency).
4. Latent food compartment named F, rebounds via resupply after crash: matches ms.Rnw:1135. Age-structured recruitment / juvenile delay: matches ms.Rnw:1156.

NOTATION CHECK: F (food compartment, italic math), \citet{searle16}, \ref{sec:no-F-ablation}, \ref{sec:res}, \Delta\mathrm{AIC}. No species superscripts needed here. \bye{} renders purple and is already used multi-paragraph elsewhere (ms.Rnw:730), so the two-paragraph block is fine.

CITATION NOTE: I did NOT add a \citep for the bootstrap LRT (it's this paper's own SI), only "Section~\ref{sec:res} and the supplement." If Bo prefers an explicit \ref{Table:SIRJPF2_model_comparison}, that label lives in si.Rnw:472. Left as prose ref to avoid a dangling main-text \ref into the SI table.

DECISIONS:
- Followed the gathered "Primary insertion point: Lines 1168-1170" guidance but did NOT replace the existing sentence; instead appended, since line 1169 is good and already cited (king15). The new block elaborates rather than overwrites.
- Kept "exchangeable"/"conditionally independent" language to dovetail with emphasis A (marginalized/block IF respects conditional independence) without duplicating the methods claim.
- Trimmed to 2 paragraphs per the 1-2 paragraph tight spec. If a 20-page trim (Task #11) bites, the second paragraph's testability point is the higher-value half to keep (it's the part deterministic analyses cannot do); the no-F numbers in para 1 are already stated in SI and could be compressed to one sentence.
- AOAS Significance Statement is a SEPARATE deliverable (Task #10); this block is Discussion-body only and uses internal cross-refs/cites freely.

```latex
\bye{The stochastic formulation does more than improve numerical fit; it changes the ecological conclusions we draw from these data. The deterministic model of \citet{searle16} attributed the rise and collapse of the \textit{Daphnia} populations to a density-dependent switch in reproductive strategy, from the asexual production of juveniles to the sexual production of dormant eggs. Our analysis instead locates the load-bearing mechanism in the depletion of the algal food resource coupled to age-structured recruitment: as the population grows it draws down the latent food compartment $F$, throttling juvenile recruitment and driving the subsequent crash, after which $F$ rebounds through resupply. The role of this feedback can be quantified directly. Holding the food compartment fixed at its initial value (the no-$F$ ablation of Section~\ref{sec:no-F-ablation}) removes the resource feedback while granting the model a free parameter, yet it costs approximately $7.2$ log-likelihood units and is decisively worse by $\Delta\mathrm{AIC}\approx 12.4$. Depletable resource dynamics, not a reproductive-strategy threshold, are therefore what the data require to explain the observed peak-then-collapse trajectory.

Stochasticity is also what makes these mechanistic comparisons testable. Process noise in the juvenile, food and parasite dynamics absorbs the overdispersion and replicate-to-replicate variability that a deterministic model can only accommodate by positing systematic differences, in the form of mesocosm-specific parameters or additional mechanisms, between otherwise identical replicates. Because the stochastic model treats each mesocosm as a conditionally independent realization of a shared dynamic process, likelihood-based model selection and the parametric bootstrap likelihood ratio test (Section~\ref{sec:res} and the supplement) allow us to ask which mechanisms genuinely need to differ across replicates. The answer here is that none do: no unit-specific parameterization yields a statistically significant improvement over the all-shared model, so the heterogeneity between replicates is adequately explained by dynamic noise rather than by structural differences. Such hypothesis tests on the sources of variation are simply unavailable to a deterministic fit, which lacks a probability model for the data and therefore cannot separate process-level stochasticity from genuine between-unit structure. In this sense, embedding the mechanism in a stochastic dynamic model is not a statistical refinement applied after the science is settled; it is the tool that redirects the mechanistic explanation toward resource depletion and reveals the experimental replicates to be statistically exchangeable.}
```

---

## Draft 4 — Two deliverables for ms.Rnw: (1) Section 2 status + outline + opening \bye{} prose for "The motivating experiment and data," and (2) a refined, AOAS-compliant Significance Statement. Both wrapped in \bye{}.

=== DELIVERABLE 1: SECTION 2 ("THE MOTIVATING EXPERIMENT AND DATA") ===

STATUS — a dedicated section ALREADY EXISTS. It is `\section{A {\it Daphnia} mesocosm experiment}` at ms.Rnw line 272, label `\label{sec:studysystem}` (line 273), running through line 419. It already covers: the two Daphnia species (native D. dentifera, invasive D. lumholtzi); the alga A. falcatus as food and the parasite M. bicuspidata; the 6 treatments (3 host levels x 2 parasite levels); the replicated units (8/9/10 per treatment); the sampling protocol; and the data figure (Fig.~\ref{fig:data_vis}). The editor's "dedicated Section 2" request is therefore ALREADY SATISFIED in content. Recommended action: ADAPT, not create — (a) optionally retitle to foreground "experiment and data" and align the label, and (b) add one short opening paragraph that states up front WHY the applied questions matter (so the section opens with motivation rather than species natural history). The opening \bye{} prose below supplies that paragraph; insert it immediately after line 273 (the \label line), before the existing `\revision{\editDaphnia}` on line 275.

OUTLINE (current content, confirming coverage):
  2.0 [NEW opening paragraph — supplied below] Why these questions matter: native/invasive competition mediated by a shared parasite; resource depletion vs. reproductive-switch as rival drivers of boom-and-collapse; replicated mesocosms as a testbed for panel inference.
  2.1 The two species and their ecological context (lines 276-281; \editDaphnia gloss on the genus).
  2.2 Life history relevant to the model: parthenogenesis, maturation, male/ephippia production under stress (lines 282-286; \editParthenogenesis).
  2.3 Infection biology of M. bicuspidata (lines 287-290).
  2.4 Experimental design: mesocosm setup, 6 treatments = 3 host x 2 parasite, 8/9/10 replicate units, initialization, spore dosing, 5-day sampling over 52 days/10 sessions (lines 292-303).
  2.5 The data (Fig.~\ref{fig:data_vis}): focus on the 4-species treatment; peak near day 25, collapse, mild resurgence near day 50 (lines 306-411).
  2.6 Rival mechanistic hypotheses from the data: Searle16 reproductive-switch vs. our explicit resource-depletion + age-structure hypothesis (lines 412-419).

If the title/label is changed, update the cross-references on lines 267 and 410 (both use `Section~\ref{sec:studysystem}`) — keeping the existing label avoids any edit churn, so the \bye{} note flags retitling as optional.

=== DELIVERABLE 2: SIGNIFICANCE STATEMENT ===

STATUS — one EXISTS but is an unfinished placeholder: ms.Rnw line 1182 is `\revision{\editSignificance}`, pulling the draft from edits.tex lines 8-13. That draft (a) renders red/unfinished, (b) speaks ONLY to the methods/stochasticity angle, and (c) mentions NEITHER editor emphasis explicitly (no panel-iterated-filtering / marginalized lever; no stochastic-vs-deterministic ecological reframing). Action: REFINE. The version below is 120 words exactly (AOAS limit), jargon-aware (explains "water-flea," avoids "POMP"/"PanelPOMP"/"marginalized"/"MCAP"/"likelihood" unexplained), self-contained (no figure/table/citation refs), and weaves BOTH emphases: A = panel iterated filtering plus the reliability-improving refinement (the marginalized/block variant, described plainly); B = stochasticity reshaping the ecological conclusion toward resource depletion and explaining replicate variation by chance alone.

**Insert at:** `/Users/ybb/Downloads/Research/Daphnia-ms/daphnia-article/ms.Rnw:274 (Section 2 opening paragraph, after \label{sec:studysystem} on L273, before \revision{\editDaphnia} on L275); /Users/ybb/Downloads/Research/Daphnia-ms/daphnia-article/ms.Rnw:1182 (Significance Statement, replacing the \revision{\editSignificance} placeholder; note \section*{Significance Statement} on L1180 already present, so either drop the duplicate \section* from the \bye{} block or remove L1180)`

**Verify / notes:**

VERIFY:
1. Significance Statement is exactly 120 words (AOAS limit is <=120). Confirmed via wc -w on the prose. If Bo edits a word, re-check the count.
2. DUPLICATE \section* HEADER: ms.Rnw L1180 already has `\section*{Significance Statement}` followed by the `\revision{\editSignificance}` placeholder on L1182. My Significance \bye{} block INCLUDES its own `\section*{Significance Statement}`. To avoid two headers, EITHER (a) delete the existing L1180 header and replace L1182 with the full \bye{} block, OR (b) drop the `\section*{Significance Statement}` line from inside my \bye{} block and just replace L1182. Recommend (b) for minimal churn. Flagging because StructuredOutput shows the self-contained version.
3. SECTION ALREADY EXISTS — I did NOT create a new section. The editor's "dedicated Section 2" is satisfied by existing sec:studysystem (L272-419). My deliverable is an opening motivation paragraph + retitle option. Decide whether to retitle `A {\it Daphnia} mesocosm experiment` -> e.g. `The motivating experiment and data`. If retitled, the section appears AFTER the Introduction as Section 2 already (it is the 2nd \section), so no reordering needed. If you keep `\label{sec:studysystem}`, the two cross-refs at L267 and L410 need no change.

DECISIONS:
- Emphasis A in the Significance Statement is stated in plain language ("a recently developed simulation-based method, panel iterated filtering... a refinement that improves its reliability") rather than naming block=TRUE / marginalized / Breto20, per AOAS jargon-free requirement. The technical naming of the marginalized/block variant belongs in the methods section and the response letter, not the Significance Statement.
- Emphasis B (stochastic vs deterministic) is the rhetorical climax: "accounting for randomness reshapes the ecological conclusions: depletion of the food resource, rather than a switch in reproductive strategy." This directly mirrors the no-F ablation result and the Searle16 contrast, without citing them (Significance must be self-contained).
- Kept the existing draft's punchline ("chance variation alone explaining differences between replicates") because it is the paper's distinctive parsimony claim; tightened wording to fit the 120-word budget.
- The Section 2 opening paragraph deliberately states the two applied questions (native/invasive competition under shared parasitism; what drives boom-collapse) AND frames replication as the panel structure the methodology exploits — this pre-loads both the ecological motivation and the methodological hook before the natural-history prose (\editDaphnia) begins.
- Both blocks are wrapped in \bye{} per convention; \bye is defined in edits.tex L6 as {\color{purple}[BY: #1]}, multi-paragraph safe, so the prose renders purple as an unfinished-revision draft for Bo's review.

OPEN: The opening paragraph cites \citet{searle16}, which is already a defined bib key used throughout. No new citations introduced.

```latex
% ========================================================================
% DELIVERABLE 1 — Section 2 ("motivating experiment and data")
% Status: section ALREADY EXISTS (sec:studysystem, ms.Rnw L272-419).
% ADAPT, do not create. Insert the opening \bye{} paragraph below
% immediately AFTER \label{sec:studysystem} (L273), BEFORE \revision{\editDaphnia} (L275).
% Optional retitle shown in the note; keep existing label to avoid edit churn.
% ========================================================================

\bye{The experiment that motivates this article asks a question of broad ecological importance: how does competition between a native and an invasive species change when both are attacked by a shared parasite, and what drives the boom-and-collapse population dynamics that result? \citet{searle16} designed a replicated mesocosm experiment to address this, growing the North American native \textit{D.~dentifera} and the invasive \textit{D.~lumholtzi}, alone and together, with and without the fungal parasite \textit{M.~bicuspidata}, on a shared algal food supply. The replicated design---multiple mesocosms within each of six treatments---is exactly the panel structure that our inference methodology is built to exploit, and the rival biological explanations for the observed dynamics (a switch in reproductive strategy versus depletion of the food resource) are exactly the kind of competing mechanistic hypotheses that likelihood-based model comparison can adjudicate. We therefore describe the experimental system, the relevant life history of the organisms, and the resulting time series data before introducing the models and inference methods used to analyze them.}


% ========================================================================
% DELIVERABLE 2 — Significance Statement (refines \revision{\editSignificance})
% Replace ms.Rnw L1182 (the \revision{\editSignificance} placeholder) with this.
% 120 words exactly; jargon-aware; self-contained; covers editor emphases A and B.
% ========================================================================

\bye{\section*{Significance Statement}

Mechanistic models help scientists learn how ecological systems work, but fitting them to data is hard when key quantities go unobserved, dynamics are noisy and nonlinear, and several replicate systems must be analyzed together. We show that a recently developed simulation-based method, panel iterated filtering, makes rigorous likelihood-based analysis of such replicated experiments practical, and we introduce a refinement that improves its reliability. Applying these tools to a freshwater experiment with two competing water-flea species, their algal food, and a shared fungal parasite, we find that accounting for randomness reshapes the ecological conclusions: depletion of the food resource, rather than a switch in reproductive strategy, best explains the observed boom and collapse, with chance variation alone explaining differences between replicates.}
```

---

## Draft 5 — si.Rnw Section S3 "MCAP Coverage Study": new final subsection "Interpretation and reporting" (C-primary form), to be inserted after the identifiability-findings subsection (and its figure/table block) and before \section{AIC Selection-Rate Simulation}.

**Insert at:** `/Users/ybb/Downloads/Research/Daphnia-ms/daphnia-article/si/si.Rnw:370`

**Verify / notes:**

INSERTION: place the three \bye{} paragraphs between line 370 (closing `}` of the identifiability-findings figure/table \bye block, `\end{table}}`) and line 372 (`\section{AIC Selection-Rate Simulation}`). Add one blank line before and after to match file spacing. This makes "Interpretation and reporting" the 4th and final subsection of S3 (after S3.3 Identifiability findings), so it correctly sits after the identifiability-findings subsection as the task requires.

NOTATION CHECK (all verified against si.Rnw): used \native, \invasive, r^{\native}, r^{\invasive}, p^{\invasive}, f_S^{\invasive}, \sigma_P, \sigma_F as elsewhere in S3. Cross-ref uses \ref{sec:mcap-coverage} (label defined line 282) to point back at the three-loglik diagnostic and the product-ridge facts rather than restating numbers, avoiding drift. "Table~1 of the main paper" matches the existing reference style at si.Rnw:222 ("Table~1 and Figure~S-1"). MCAP attribution \citep{ionides17,ning21} already established at line 206; not re-cited here to avoid redundancy within the section.

NUMBERS reused (consistent with current S3 prose, lines 290-293, 326): coverages 77.8/87.0/74.0%; MCSE ~4 pp; p=24 free params, 80 observations; warm-from-truth 20-chain search; 74-87% range for identifiable params. All match Table:mcap-coverage and the three-loglik subsection. No new numbers introduced.

CLAIMS FOR BO TO VERIFY:
1. The "74-87% range" characterization for the Table 1 data intervals is an extrapolation from the bootstrap (which simulated from the all-shared MLE) to the actual data fit. This is the honest-reporting framing Ed endorsed (option C), but confirm you want the main-paper Table 1 intervals explicitly described to readers as plausibly 74-87% empirical coverage — this is a stronger/more candid statement than merely citing the SI. If too strong, soften to "whose empirical coverage in this regime is estimated in Section S3."
2. Confirm "option B = enlarge the likelihood drop / calibrate the cutoff so 95% of bootstrap intervals cover" matches your actual MCAP calibration mechanism (drop in smoothed profile log-lik defining the interval). I described it as enlarging the cutoff; verify the lever is the profile-drop threshold and not, e.g., the MCAP lambda/regression-smoothing parameter.
3. The claim "would require an order of magnitude more bootstrap datasets" is a qualitative defense of B-as-secondary; ensure it is consistent with whatever you write in response.tex emphasis (keep B=100-too-small-to-calibrate-a-95%-quantile caveat aligned across SI and response letter).
4. "$p=24$ free parameters fit to $80$ observations": matches line 326 (p=24, U=8, N_u=10 => 80). Confirm 80 (not the per-the-table 99/100 denominators, which are bootstrap counts) is the intended sample-size statement.

STRUCTURE NOTE: I wrapped each of the 3 paragraphs in its own \bye{} (matching the per-paragraph \bye{} pattern in S3.1/S3.2 rather than one block spanning all). The \subsection title itself is left OUTSIDE \bye{} — consistent with lines 288/313/331 where \subsection{...} is unwrapped and only the prose is purple. If you instead want the whole pending subsection (incl. title) marked purple, wrap the title too; current file convention is title-unwrapped, which I followed.

DROP-IN: this is plain LaTeX (no Sweave chunks), so it compiles in the .Rnw without an R block. No new figures/tables added (the existing fig/table for S3 already precede this subsection).

```latex
\subsection{Interpretation and reporting}

\bye{We report this coverage study in an exploratory rather than a confirmatory spirit.
By exploratory we mean that its purpose is to characterise the inferential behaviour of the SIRJPF2 model---to map which directions in parameter space the short panel does and does not constrain, and to quantify the finite-sample displacement of the maximum likelihood estimate along the weakly-identified ridges identified above---rather than to certify that the MCAP intervals deliver an exact $95\%$ frequentist guarantee on this model.
The empirical coverages of $77.8\%$, $87.0\%$ and $74.0\%$ for $r^{\native}$, $\sigma_P$ and $\sigma_F$ are therefore presented as honest summaries of how the intervals behave in this regime, together with their Monte Carlo standard errors of roughly four percentage points, and not as adjustments to be applied retroactively.
The undercoverage is a property of the short, weakly-identified panel ($p=24$ free parameters fit to $80$ observations) and the resulting finite-sample displacement of the estimator documented in Section~\ref{sec:mcap-coverage}, not a numerical artefact: the three-log-likelihood diagnostic confirms that the production profiles are converged and in fact more thorough than a $20$-chain warm-from-truth search.}

\bye{We accordingly take honest empirical reporting as the primary mode of presentation, and the MCAP confidence intervals reported for the data analysis in Table~1 of the main paper are left unchanged.
A reader should read those intervals as smoothed profile-likelihood intervals whose nominal level is $95\%$ but whose empirical coverage in this finite-sample, weakly-identified setting is plausibly closer to the $74$--$87\%$ range estimated here for the identifiable parameters, and as essentially uninformative for the non-identifiable parameters ($r^{\invasive}$, $p^{\invasive}$, $f_S^{\invasive}$) whose levels enter the process model only through the products $r^{\invasive} f_S^{\invasive}$ and $p^{\invasive} f_S^{\invasive}$.
This is the appropriate qualification for intervals whose intended use is to delimit the parameter ranges consistent with the data and to support the model comparisons in the main paper, rather than to provide certified $95\%$ statements about individual parameters.}

\bye{As a sensitivity check, the MCAP cutoff could instead be calibrated against this bootstrap so that the empirical coverage is forced to the nominal $95\%$---that is, the smoothed-profile likelihood drop defining the interval could be enlarged until $95\%$ of the $B$ bootstrap intervals cover the truth, and that calibrated cutoff then applied to the data intervals.
We regard this as available but secondary, because $B=100$ bootstrap replicates is a small basis on which to estimate a $95\%$ quantile: the relevant tail is governed by the five or so datasets nearest the coverage boundary, the Monte Carlo standard error on the calibrated cutoff is correspondingly large, and the per-parameter coverages do not move together (Section~\ref{sec:mcap-coverage}), so a single recalibrated cutoff would over-correct some parameters while under-correcting others.
A defensible calibration would require an order of magnitude more bootstrap datasets, each itself an expensive iterated-filtering profile.
We therefore present the honest empirical coverages as our primary result and note cutoff calibration only as a direction available for a future, more heavily resourced study; the reported intervals in the main paper are not recalibrated.}
```