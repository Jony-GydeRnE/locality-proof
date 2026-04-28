# n=8 singleton Laurent cascade

> **Headline: 100 / 100 step-1 survivors at $n=8$ die via a depth-1
> Laurent cascade. Verified symbolically with sympy. Zero failures.**

## What the test did

For each of the 100 non-triangulation step-1 survivors $M$ at $n=8$,
the script searches for a triple $(Z, k, \mathrm{fp})$ — a kill zone
$\mathcal Z_{r, r+2}$, a Laurent order $k$, and a free-variable
fingerprint $\mathrm{fp}$ — such that:

1. **$M$ contributes** to the order-$X_S^{-k}$ Laurent expansion of
   $1/\prod_{c \in M} X_c$ on zone $Z$ with nonzero scalar on
   $\mathrm{fp}$.
2. **Every other multiset** appearing in the same fingerprint equation
   is step-1 killable at *some* zone (not necessarily $Z$).

If both hold, substituting the cousins' zeros into the fingerprint
equation collapses it to $\mathrm{scalar}(M) \cdot a_M = 0$, hence
$a_M = 0$.

The candidate generation mimics the worked $n=7$ recipe: for each
substitute chord $U \in M$ at zone $Z$, take its companion $Y_U$ as
the fingerprint numerator and the rest of $M$'s free chord factors
as the denominator; try Laurent orders $k = \mathrm{leading} + 1$
and $k = \mathrm{leading} + 2$.

## How long it took

- Smoke test (5 survivors): 5.5 min wall-clock
- **Full run (100 survivors): ~76 min wall-clock (4 536 s)**
- Per-orbit audit (13 orbits, all (Z, U) recipes per orbit): ~10 min

## Commands

```bash
cd scripts
python3 cascade_kill_n8.py            # full singleton-cascade run (~76 min)
python3 cascade_kill_n8.py 5          # smoke test on first 5 survivors
python3 analyze_recipes.py            # per-survivor analysis (seconds)
python3 audit_orbits.py               # per-orbit (Z, U) audit (~10 min)
```

## What the result means

### Singleton cascade kills (the main verifier)

| Laurent order $k$ | # survivors |
|---|---:|
| 3 (= leading + 1) | 90 |
| 4 (= leading + 2) | 10 |

All 100 die at depth 1. Every cousin in every fingerprint equation
is step-1 killable at some neighbouring zone — exactly the cascade
pattern proven at $n = 7$ for the seven fish.

### Cyclic orbit decomposition (`outputs/orbits.md`)

The 100 survivors break into **13 cyclic orbits**: 12 of size 8
(action of $\mathbb Z_8$ is free) and 1 of size 4 (stabilised by the
$\mathbb Z_2$ subgroup).

### Per-orbit audit (`outputs/orbit_audit_summary.md`)

For each orbit, every K2-clean $(Z, U)$ pair was tested. Result:

> **Every orbit has at least one valid depth-1 $(Z, U)$ recipe; total
> 26 valid recipes across 13 orbits.**

Many orbits have multiple valid recipes (orbits 6 and 12: four each),
which gives slack for any future combinatorial proof of universality.

### Φ-I rule status (`findings.md`)

Φ-I (the "frame" rule: kill at zone $\mathcal Z_{r, r+2}$ when both
endpoints $\{r, r+2\} \subseteq V_\text{missed}(M)$) matches 48/100
survivors. **Φ-I is therefore not the universal rule.** The
audit-derived rule that *every* K2-clean $(Z, U)$ pair with
$\ell_Z \ge 2$ produces a recipe whose cousins all step-1-die does
match all 100 — but proving that universally is the open piece.

## Files

```
n8/
├── README.md                              ← this file
├── findings.md                            ← per-orbit headline takeaways
├── scripts/
│   ├── cascade_kill_n8.py                 ← main singleton-cascade verifier
│   ├── analyze_recipes.py                 ← per-survivor analysis
│   └── audit_orbits.py                    ← per-orbit (Z, U) enumeration
└── outputs/
    ├── results_cascade_n8.txt             ← full cascade trace (100 blocks)
    ├── recipe_analysis.md                 ← per-survivor markdown table
    ├── orbits.md                          ← orbit decomposition
    ├── analysis_summary.txt               ← short summary
    ├── orbit_audit.json                   ← per-orbit audit, machine-readable
    └── orbit_audit_summary.md             ← per-orbit audit, human-readable
```

## What remains unknown at n=8

- A *combinatorial proof* (not just an empirical check) of the audit-
  derived rule "every K2-clean $(Z, U)$ with $\ell_Z \ge 2$ closes".
- The asymptotic rate of survivor growth as $n$ increases (data side
  of Item 4; lives in `paper/step-one-kill-statistics/`).

The big *known unknown* — the n=8 → n=9 transition where singleton
cascade ceases to be universal — is documented in the sibling [`n9/`](../n9/)
folder.
