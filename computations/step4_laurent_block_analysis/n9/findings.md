# n=9 findings — quick reference

For the full block-rule kill writeup see [`cluster_findings.md`](cluster_findings.md);
for per-anomaly verdicts see [`anomaly_summary.md`](anomaly_summary.md);
for per-failure detail see the [`diagnostics/`](diagnostics/) folder.

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

> **Block-rule Laurent kill works at $n=9$.**

The singleton cascade rule that closed every $n=7, 8$ survivor
sometimes fails at $n=9$, but the depth-1 fingerprint matrix on
the cluster around orbit 22 has trivial cluster nullspace, so the
4 failures (and their cluster) are killed jointly.

## What changes for the proof

The all-$n$ conjecture (Item 10 of the master plan) was originally
stated as *"every step-1 survivor admits its own depth-1 cascade
recipe."* The $n=9$ data refines this to:

> *the depth-1 fingerprint matrix on each connected component of
> the depth-1 cousin graph has trivial cluster nullspace.*

The singleton rule is a special case (cluster of size 1).
