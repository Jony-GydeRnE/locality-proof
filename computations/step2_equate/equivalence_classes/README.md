# Step-2 equivalence-class analysis

**Script:** `step2_classes.py` (parametric in n).

## Question
At each cyclic zone Z_{r,r+2}, the bare relation X_{r+1,r-1} = −X_{r,r+2} forces equalities of the form a_M = a_{M'} between coefficient multisets that differ by one bare↔special swap. Cyclic shifts give more equalities. Together they partition the survivors into **equivalence classes**.

If all C_{n−2} triangulations land in ONE class, that's the unitarity statement.

## Two modes (both reported)
- **PERMISSIVE** — every bare↔special swap (ignoring multi-term partial-fraction effects). Over-equates and incorrectly forces triangulations to 0 via false-positive links to Step-1-killed siblings. Reported as a *diagnostic*.
- **CONSERVATIVE** — only swap edges where BOTH endpoints are Step-1 survivors. A *lower bound* on the true Step-2 equivalence (any indirect chain through a non-survivor intermediate is suppressed).

## Headline numbers (conservative mode)
| n | tri survivors | non-tri survivors | classes (cons.) | cyclic orbit-shapes |
|---|---:|---:|---:|---:|
| 5 |  5 |  0 |  1 (size 5) | 1 orbit (the unitarity class) |
| 6 | 14 |  0 |  5 (sizes 1, 1, 4, 4, 4) | 2 orbits (sizes 2 + 3) |
| 7 | 42 |  7 | 21 (7×size-1, 7×size-2, 7×size-4) | 3 orbits, each of size 7 |

So conservative gives full unitarity at n=5 but fragments at n≥6. The missing equivalences (which the existing `full_ansatz_nullspace/` scripts confirm via nullspace dim = 1) come from indirect chains through non-survivor intermediates that the conservative rule blocks.

## Run
```bash
python3 step2_classes.py 7
```
