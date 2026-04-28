# §7 — Step-3: block-rule cascade kill at $n=9$

**Backs paper §7. Closes locality at $n=9$ from the cyclic 1-zeros
alone.**

## What this paper proves

> **Locality at $n=9$.** If a rational ansatz $B$ of mass dimension
> $-2(n{-}3)$ with at most simple planar poles satisfies the nine
> cyclic hidden-zero conditions $B|_{\mathcal Z_r}=0$, then $a_M=0$
> for every non-triangulation size-$6$ multiset $M$.

This is unconditional — the non-local half of the
$\mathrm{Tr}(\phi^3)$ uniqueness theorem at $n=9$, proved end-to-end
from the 1-zero conditions.

Combined with the $d$-subset / homogeneity-rescaling unitarity
argument from `paper/Details missing in Hidden zero -> unitarity
Rodina proof/` and notes #11, #18, this gives the full
$B=c\cdot A_9^{\text{tree}}$ uniqueness statement at $n=9$.

## The mechanism

| $n$ | Step-1 non-tri survivors | Cluster picture | Layer needed |
|---:|---:|---|---|
| 4, 5, 6 | 0 | none — Step-1 closes locality | Step-1 |
| 7 | 7 (1 orbit) | 1 singleton cluster | Step-1 + singleton cascade |
| 8 | 100 (13 orbits) | 13 singleton clusters | Step-1 + singleton cascade |
| 9 | 1011 (113 orbits) | 23 singletons + **1 non-trivial 90-orbit cluster** | Step-1 + singleton cascade + **block-rule** |

At $n \le 8$ every cluster in the depth-1 cousin graph is a
singleton: each orbit's depth-1 fingerprint equation has all
non-local cousins step-1-killable, and the equation collapses to
$a_M=0$ on its own.

At $n=9$ singleton cascade *sometimes* fails: 4 orbit
representatives (22, 46, 88, 108 of 113) have no depth-1 recipe
whose every cousin is step-1-killable, because their cousins are
themselves Step-1 survivors in *other* orbits. This is the
**coupled subsystem** phenomenon, and it requires the **block-rule
kill**: the depth-1 fingerprint matrix on the entire connected
cluster has trivial cluster nullspace, forcing every cluster orbit
to vanish jointly.

## n=9 numerical anchor

Verified at $n=9$ on the 90-orbit cluster around orbit 22:

- 506 depth-1 fingerprint equations
- 90 cluster columns + 53 external columns (32 triangulations + 21
  non-cluster non-tri survivors)
- **Cluster rank = 90 = full**; cluster residual $r = 0$
- Setting non-cluster non-local survivor columns to 0 (each killed
  by its own single-orbit cascade) and any consistent triangulation
  assignment, the cluster columns are forced to 0.

The 23 non-cluster orbits are each killed by their own single-orbit
depth-1 cascades (verified, including a re-run of the originally
timed-out orbit 68 which closes in 224 s with all 7 cousins
step-1-killable). Combined with Step-1 covering 904,752 / 905,763
non-triangulation multisets directly, every non-local coefficient
at $n=9$ vanishes.

## Files in this folder

- `step3_block_rule_n9.tex` — the proof at $n=9$.
- `step3_block_rule_n9.pdf` — compiled.

## Where the supporting computations live

- **Consolidated per-orbit locality audit** (one-stop reference):
  `../../computations/step4_laurent_block_analysis/n9/outputs/n9_locality_status.md`
  This file accounts line-by-line for each of the 113 orbit
  representatives, marking it as cluster or non-cluster and showing
  its kill mechanism (block-rule or single-orbit cascade recipe).
- Cluster matrix construction at $n=9$:
  `../../computations/step4_laurent_block_analysis/n9/`
- Single-orbit cascade results (per orbit):
  `../../computations/step4_laurent_block_analysis/n9/outputs/results_cascade_n9_reps.json`
- Per-orbit anomaly diagnostics:
  `../../computations/step4_laurent_block_analysis/n9/diagnostics/`
- Step-2 flip-graph at $n=9$ (16 components on 49 triangulation orbits):
  `../../computations/step2_equate/flip_graph_n9/`
  (This is a *unitarity* artifact; not part of the locality proof
  here. See "Locality vs unitarity" below.)

## Locality vs unitarity: a clarification

The `flip_graph_n9` analysis identified 4 of the 16 Step-2
triangulation components as not bridged by the cluster matrix's
external columns (the "untouched components" with chord-length
signature (2,2,3,3,4,4) — fan-class triangulations).

**These 4 components are a unitarity question, not a locality
question.** Triangulations are local terms — they should not be
killed; they are designed to survive Step-1 (see the local-survival
lemma in `paper/step1 Kill technique and statistics/step-one-doesnt-kill-triangulations/`).
What the unitarity argument needs is to equate all 49 triangulation
$c$-values to a single common scalar. That is handled by the
$d$-subset argument in `paper/Details missing in Hidden zero ->
unitarity Rodina proof/details for Rodina D subset argument/` and
notes #11, #18 — independently of the cascade machinery here.

## What this conjectures for general $n$

For every $n$, every cluster $\mathcal C$ of Step-1 survivor orbits
has cluster-column rank $\ge |\mathcal C|-1$ on its depth-1
fingerprint matrix; if a residual $r\ge 1$ ever appears, it is
closed by $r$ higher-Laurent-order anchor equations. Singleton
cascade is the special case $|\mathcal C|=1$, $r=0$; the block-rule
generalises it to clusters of any size. Verified at $n=4$–$9$:
$r=0$ in every observed case.

## What remains for the all-$n$ locality theorem

1. **Cluster finiteness** — every cluster is finite at every $n$.
   Likely follows from finiteness of size-$(n{-}3)$ chord
   multisets; record as a lemma.
2. **Cluster rank lower bound** —
   $\mathrm{rank}(M_{\mathrm{FP}}^{\mathcal C})\ge |\mathcal C|-1$
   on cluster columns for all $n$. This is the core
   combinatorial-geometric content.
3. **Anchor existence when $r=1$** — through $n=9$ no cluster has
   shown $r\ge 1$, but the all-$n$ proof should accommodate the
   possibility (an explicit higher-Laurent anchor closing the
   residual).

These are purely locality obligations. Unitarity is independent
and handled by the $d$-subset argument.
