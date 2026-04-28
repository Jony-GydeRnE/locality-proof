# n=9 locality — fully proven

> **Locality at $n=9$ is fully proven.** Every non-triangulation
> coefficient is forced to zero by Step-1 + 23 single-orbit depth-1
> cascades + a block-rule kill on a 90-orbit coupled cluster.

## What turned out to be happening

What looked like four stubborn "non-local survivors" at $n=9$ — orbits
22, 46, 88, 108, with no private depth-1 cascade recipe — turned out to
be **interrelated**. They sit inside a 90-orbit coupled cluster of the
depth-1 cousin graph and **get killed together**: the cluster's depth-1
fingerprint matrix has rank 90 (full) and nullity 0, so every cluster
coefficient is forced to zero jointly.

## Headline numbers

| | Value |
|---|---:|
| Total size-6 chord multisets | 906 192 |
| Triangulations ($C_7$, local; never killed) | 429 |
| Non-triangulations | 905 763 |
| **Killed by Step-1 alone** | **904 752** (99.889% of non-tri) |
| Step-1 non-tri survivors | 1 011 → 113 cyclic orbits |
| Killed by single-orbit depth-1 cascade | 23 orbits |
| Killed jointly by block-rule cluster (rank 90) | 90 orbits |

## What to read

- [`findings.md`](findings.md) — one-page summary of the proof.
- [`outputs/n9_locality_status.md`](outputs/n9_locality_status.md) — **per-orbit consolidated audit (all 113 reps).** This is the file the paper cites.
- [`outputs/cluster_partial.md`](outputs/cluster_partial.md) / [`cluster_partial.json`](outputs/cluster_partial.json) — block-rule cluster matrix (90 cluster columns, 506 rows, rank 90, nullity 0).
- [`outputs/orbits_n9.md`](outputs/orbits_n9.md) — orbit manifest (113 reps).
- [`outputs/results_cascade_n9_reps.txt`](outputs/results_cascade_n9_reps.txt) — per-rep singleton-cascade trace.
- [`diagnostics/`](diagnostics/) — per-orbit deep-search reports cited from the paper (orbits 22, 25, 28, 34, 35, 46, 68, 88, 108).

## Pipeline (run in order)

```bash
cd scripts
python3 orbit_decomposition.py            # 1 011 step-1 survivors → 113 cyclic orbits
python3 cascade_kill_orbit_reps_n9.py     # singleton depth-1 cascades, one per orbit rep
python3 cluster_analysis.py 22            # BFS-grow the depth-1 cousin cluster from orbit 22
python3 cluster_partial_analysis.py       # rank/nullity of the 90-orbit cluster fingerprint matrix
python3 build_locality_status.py          # consolidated per-orbit audit → outputs/n9_locality_status.md
```

## Pre-block-rule diagnostics (retired)

Earlier scratch work — including a misleadingly named
`n9_failure_equations_readable.txt` whose "failures" are exactly the
4 orbits the block rule was invented to handle jointly — lives at
[`old stuff/n9 pre-block-rule diagnostics/`](../../../old%20stuff/n9%20pre-block-rule%20diagnostics/).
