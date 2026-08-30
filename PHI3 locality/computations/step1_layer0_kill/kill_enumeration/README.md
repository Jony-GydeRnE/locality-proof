# Step-1 (Layer-0) kill enumeration

**Script:** `step1_uncaught.py` (parametric in n).

## Question
At each of the n cyclic zones Z_{r,r+2}, the leading-order Laurent argument (special chord X_{r,r+2} → ∞) kills certain coefficient multisets. Across all n zones, **which multisets remain uncaught?**

## Output
- Console report: total multisets, # uncaught, # of those that are triangulations vs non-triangulations.
- PDF in `outputs/`: every non-triangulation survivor drawn as a filled n-gon with red chords.

## Headline numbers
| n | total multisets | uncaught (= tri + non-tri) | non-tri survivors |
|---|---:|---:|---:|
| 5 | 15 | 5 + 0 | **0** |
| 6 | 165 | 14 + 0 | **0** |
| 7 | 2 380 | 42 + 7 | **7** |
| 8 | 42 504 | 132 + 100 | **100** |
| 9 | 906 192 | 429 + 1 011 | **1 011** |

The triangulation counts are exactly the Catalan numbers C_{n−2} as expected (the tree amplitude survives every zone). The non-triangulation survivors are the "diagrams that refuse to die" — the inputs to the Step-2 equality-chain analysis.

## Run
```bash
python3 step1_uncaught.py        # prompts for n
python3 step1_uncaught.py 8      # n on the CLI
```
