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

### Visualization galleries (the 90 cluster orbits, three views)

These let you see the 90 cluster unknowns and the structure of the 506 equations directly, instead of only as numbers.

- [`outputs/cluster_90_gallery.pdf`](outputs/cluster_90_gallery.pdf) — one diagram per cluster orbit (orbit IDs 1..113, filtered to the 90). Black 9-gon + red chord overlays. Each diagram is one column of the rank-90 cluster system. The 4 cascade-failure seeds (orbits 22, 46, 88, 108) are flagged `[SEED]`. Built by [`scripts/cluster_gallery.py`](scripts/cluster_gallery.py).
- [`outputs/cluster_dihedral_gallery.pdf`](outputs/cluster_dihedral_gallery.pdf) + [`outputs/dihedral_reduction.md`](outputs/dihedral_reduction.md) / [`.json`](outputs/dihedral_reduction.json) — the same 90 orbits modded out by reflection. **90 Z_9 orbits collapse to 55 D_9 groups**: 7 palindromic (reflection-self-symmetric, alone in their group), 35 reflection pairs (two Z_9 orbits merged), 13 "palindromic-leak" (alone in cluster but reflection partner sits OUTSIDE the cluster — the cluster is *not* closed under D_9). Sanity: 7 + 13 + 2·35 = 90. Built by [`scripts/dihedral_reduction.py`](scripts/dihedral_reduction.py).
- [`outputs/row_structure.pdf`](outputs/row_structure.pdf) / [`row_structure.json`](outputs/row_structure.json) — visualizes the 506 equation rows: per-orbit row counts, per-(orbit, zone) heatmap (sorted both by orbit ID and by D_9 group), and a count-level reflection-symmetry test. **Finding:** even at the count level the system is *not* fully D_9-invariant (only 6/35 reflection pairs and 2/7 palindromes have matching per-zone equation counts). The kill mechanism privileges the +1 vertex direction at each zone, so the cluster matrix has Z_9 symmetry but breaks D_9. Built by [`scripts/row_structure_viz.py`](scripts/row_structure_viz.py).

The visualization scripts read only `cluster_partial.json` + `orbits_n9.json` and don't require the 10-hour build to be redone.

### Saving the matrix entries (for future block-heatmap viz)

`cluster_partial_analysis.py` now also writes [`outputs/cluster_partial_matrix.npz`](outputs/cluster_partial_matrix.npz) with the integer entries of the 506×143 fingerprint matrix. The previous run's JSON dropped these to keep file size small; re-running the script regenerates the matrix once the .npz is needed for block-coupling visualization.

## Pipeline (run in order)

```bash
cd scripts
python3 orbit_decomposition.py            # 1 011 step-1 survivors → 113 cyclic orbits
python3 cascade_kill_orbit_reps_n9.py     # singleton depth-1 cascades, one per orbit rep
python3 cluster_analysis.py 22            # BFS-grow the depth-1 cousin cluster from orbit 22
python3 cluster_partial_analysis.py       # rank/nullity of the 90-orbit cluster fingerprint matrix (slow; ~10h)
python3 build_locality_status.py          # consolidated per-orbit audit → outputs/n9_locality_status.md

# Visualization (fast; reads only the saved JSONs):
python3 cluster_gallery.py                # 15-page PDF, one diagram per cluster orbit (90 orbits)
python3 dihedral_reduction.py             # 90 Z_9 orbits → 55 D_9 groups + mini-gallery
python3 row_structure_viz.py              # row counts, zone heatmap, D_9 symmetry test of equations
```

## Pre-block-rule diagnostics (retired)

Earlier scratch work — including a misleadingly named
`n9_failure_equations_readable.txt` whose "failures" are exactly the
4 orbits the block rule was invented to handle jointly — lives at
[`old stuff/n9 pre-block-rule diagnostics/`](../../../old%20stuff/n9%20pre-block-rule%20diagnostics/).
