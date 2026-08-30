# n=9 findings — quick reference

For the per-orbit consolidated audit see [`outputs/n9_locality_status.md`](outputs/n9_locality_status.md);
for per-orbit deep-search reports see the [`diagnostics/`](diagnostics/) folder.

This file is a one-page summary.

## Numbers

| | Value |
|---|---|
| Step-1 survivors at $n=9$ | 1 011 |
| Cyclic orbits | 113 (112 size-9 + 1 size-3) |
| Singleton-cascade reps that succeeded | 94 |
| Reps that timed out | 15 |
| Reps that genuinely failed | **4** = {22, 46, 88, 108} |
| BFS cluster from orbit 22 | 90 orbits |
| Cluster $M_{FP}$ rank / nullity | **90 / 0** |
| Cluster $M_{FP}$ residual $r$ | **0** |

## Verdict

> **Locality at $n=9$ is FULLY PROVEN.**
>
> Every non-triangulation coefficient is forced to zero by Step-1 +
> the 23 non-cluster single-orbit cascades + the block-rule kill on
> the 90-orbit cluster.

The singleton cascade rule that closed every $n=7, 8$ survivor
sometimes fails at $n=9$ (orbits 22, 46, 88, 108), but the depth-1
fingerprint matrix on the 90-orbit cluster has trivial cluster
nullspace, so the 4 failures and their cluster are killed jointly.

UNITARITY (all triangulation $c$-values equal a common scalar) is a
**separate** question, handled by the d-subset argument in
`../../paper/Details missing in Hidden zero -> unitarity Rodina proof/details for Rodina D subset argument/`
and `notes/18. d subset rigorous proof.pdf`.

## What changes for the locality conjecture

The all-$n$ conjecture (Item 10 of the master plan) was originally
stated as *"every step-1 survivor admits its own depth-1 cascade
recipe."* The $n=9$ data refines this to:

> *the depth-1 fingerprint matrix on each connected component of
> the depth-1 cousin graph has trivial cluster nullspace.*

The singleton rule is a special case (cluster of size 1).
