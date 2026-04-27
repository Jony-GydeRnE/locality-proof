# Dual experiment — "X₁₃ never special" (Item 6)

**Script:** `step1_dual.py` (parametric in n; default excludes zone Z₁,₃).

## Question
What if we run the Step-1 kill across **all zones except one**? Which multisets become "extra survivors" — the ones that the excluded zone uniquely killed?

This isolates the contribution of one specific zone to the overall Step-1 kill, and identifies the multisets that the Step-2 equality-chain mechanism must handle when that zone is taken off the table.

## Output
- Console: baseline survivor count (all zones) vs dual survivor count (n−1 zones), with the EXTRAS = uniquely-killed-by-the-excluded-zone multisets.
- PDF in `outputs/`: the extras drawn as n-gons with red chord multisets.

## Headline numbers (excluding Z₁,₃, so X₁₃ never special)
| n | baseline survivors | dual survivors | EXTRAS (= uniquely killed by Z₁,₃) |
|---|---:|---:|---:|
| 7 | 49  (= 42 tri + 7 non-tri)   | 70    (= 42 tri + 28 non-tri)  | **21** (all non-tri) |
| 8 | 232 (= 132 tri + 100 non-tri)| 404   (= 132 tri + 272 non-tri)| **172** (all non-tri) |
| 9 | 1 440 (= 429 tri + 1 011 non-tri) | 2 828 (= 429 tri + 2 399 non-tri) | **1 388** (all non-tri) |

Two sanity facts the experiment confirms:
1. The Catalan triangulation count is invariant to excluding a zone.
2. All extras are non-triangulations (no tree-amplitude diagram is ever uniquely killed by a single zone — every zone's kills can be covered by the other n−1 via cyclic symmetry, *but only on the triangulation set*).

## Run
```bash
python3 step1_dual.py 8                # default excludes Z_{1,3}
python3 step1_dual.py 8 --exclude 3    # excludes Z_{3,5} instead
```
