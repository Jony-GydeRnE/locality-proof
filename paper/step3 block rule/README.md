# §7 — Step-3: block-rule cascade kill

**Backs paper §7. New section for the all-$n$ direction.**

## The setup

At $n = 7, 8$ the Laurent cascade kill is **singleton**: each step-1
survivor has its own depth-1 cascade recipe whose every cousin is
step-1 killable at one neighbouring zone, and the equation collapses to
$a_M = 0$.

At $n = 9$ singleton cascade *sometimes* fails. Four orbit
representatives (22, 46, 88, 108 of 113) have no depth-1 recipe whose
every cousin is step-1 killable — because their cousins are themselves
step-1 survivors in *other* orbits. This is the **coupled subsystem**
phenomenon.

## What this file proves

The block-rule kill: **the depth-1 fingerprint matrix on each
connected component of the depth-1 cousin graph has trivial cluster
nullspace.** Equivalently, taking many depth-1 fingerprint equations
together — not one at a time — collectively forces every cluster
orbit's coefficient to vanish.

Verified at $n = 9$ on the 90-orbit cluster around orbit 22:
- 506 depth-1 fingerprint equations
- 90 cluster columns + 53 external columns (32 triangulations + 21
  non-cluster non-tri survivors)
- **Cluster rank = 90 = full**, cluster residual $r = 0$
- Setting external columns to 0 forces every cluster orbit to vanish,
  *including* the 4 "failure" orbits killed jointly with the cluster.

## Files

- `step3_block_rule_n9.tex` — the proof at $n = 9$.
- `step3_block_rule_n9.pdf` — compiled.

## Where the supporting computations live

- Cluster matrix construction at $n = 9$:
  `../../computations/step4_laurent_block_analysis/n9/`
- Per-orbit anomaly diagnostic:
  `../../computations/step4_laurent_block_analysis/n9/scripts/diagnose_anomalies.py`
- Step-2 flip-graph connectivity at $n = 9$:
  `../../computations/step2_equate/flip_graph_n9/`
  (16 components, 12 touched by cluster, 4 untouched — the remaining
  unitarity gap).

## What remains open at n=9

The block-rule kill closes the 90-orbit cluster, but the Step-2
flip-graph on the 49 triangulation orbits has 16 components. Of those,
12 are bridged by cluster external columns; **4 are not.** Bridging
these 4 is the central open question for $n = 9$ unitarity, and the
file flags it explicitly.

## What this conjectures for general n

For every $n$, the depth-1 fingerprint matrix on each connected
component of the depth-1 cousin graph has trivial cluster nullspace.
Singleton cascade is the special case of cluster size 1; block kill
is the generalisation needed when clusters grow.
