# n=9 Step-2 flip-graph data

This folder contains the Step-2 bare/special swap data for the
triangulations of the 9-gon, plus exploratory scripts that look at how
the flip-graph components relate to the n=9 cluster matrix.

## Important framing

> **Locality at $n=9$ is fully proven** by the cluster matrix
> (rank 90 = full on the 90-orbit cluster) plus the 23 non-cluster
> single-orbit cascades plus Step-1 directly — every non-triangulation
> coefficient at $n=9$ is forced to zero. See
> [`../../step4_laurent_block_analysis/n9/`](../../step4_laurent_block_analysis/n9/).
>
> **Triangulations are LOCAL terms.** Step-2 swap data is about
> *equating* triangulation coefficients (a unitarity question), not
> about killing them. The 4 "untouched" Step-2 components surfaced
> by the cluster cross-check were a misframing — they were never a
> locality gap. Whether all triangulation $c$-values reduce to a single
> common scalar at $n=9$ is handled by a separate **d-subset argument**
> (see `notes/18.` and
> `../../paper/Details missing in Hidden zero -> unitarity Rodina proof/`),
> independent of the cascade machinery.

## What's in this folder

| File | What it is |
|---|---|
| `flip_graph_n9.py` | Verifier: enumerates 429 triangulations, computes cyclic-orbit decomposition (49 orbits), Step-2 swap edges (378 at triangulation level / 41 at orbit level), and connected components (16). |
| `outputs/n9_triangulation_orbits.json` | 49 orbits + canonical reps + sizes. |
| `outputs/n9_step2_bare_swap_pairs.json` | every Step-2 swap edge at triangulation level + orbit level + which zones produce it. |
| `outputs/n9_flip_graph_connectivity.txt` | human-readable verdict, including per-component representatives. |
| `outputs/untouched_components_readable.md` | structural data on the 4 components ({1, 3, 7, 14}) that the cluster matrix's external triangulation columns don't reach. *Used for descriptive characterization only — not a locality gap.* |
| `diagrams/component_<id>_orbit_<id>.png` | PNG 9-gon diagrams for each orbit in those 4 components. |
| `visualize_untouched.py` | script that produces the readable + diagram outputs. |
| `build_bridges.py`, `find_bridge_equations.py` | exploratory scripts that searched for triangulation-bridging equations across Step-2 components. **NOT a locality test.** Kept for reference; results superseded by the d-subset argument cited above. |

## Headline numbers

| | Value |
|---|---:|
| Triangulations of the 9-gon | 429 = $C_7$ ✓ |
| Cyclic orbits | 49 (47 of size 9, 2 of size 3) |
| Step-2 swap edges (triangulation level) | 378 |
| Step-2 swap edges (orbit level) | 41 |
| Connected components (orbit level) | 16 |
| Components touched by cluster matrix's external triangulation columns | 12 |
| Components NOT touched ("fan-class", chord-length signature (2,2,3,3,4,4)) | 4 → orbits {1, 2, 27, 45}, {4, 5, 11, 49}, {12, 13, 32, 48}, {29, 30, 33, 47} |

## Why we walked through this exploration

The cluster cross-check produced a tantalising 12-of-16 pattern, and we
chased the question of whether depth-1 / depth-2 fingerprint equations
could bridge the remaining 4 components into the cluster's reach. The
honest answer (after re-reading the goal) is that **bridging
triangulation components is a unitarity question, not a locality one**,
and unitarity is established separately. So this folder records the
flip-graph data and the structural observation that the 4 untouched
components share signature `(2,2,3,3,4,4)` — interesting combinatorics
on its own, but not a locality gap to be closed here.

## Cross-references

- Locality proof at $n=9$ (cluster matrix + cascades):
  `../../step4_laurent_block_analysis/n9/findings.md` and
  `../../step4_laurent_block_analysis/n9/outputs/n9_locality_status.md`
- Triangulation equating (unitarity via d-subset):
  `../../../paper/Details missing in Hidden zero -> unitarity Rodina proof/details for Rodina D subset argument/`
  and `../../../notes/18. d subset rigorous proof.pdf`
- Step-2 conservative classes at $n=7$ (sister folder):
  `../equivalence_classes/`

## Run

```bash
python3 flip_graph_n9.py            # the only test that's "load-bearing" here
python3 visualize_untouched.py      # generates readable + PNG diagrams
```
