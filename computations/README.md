# `computations/` — experiments backing each proof step

This folder is organized by the **logical steps of the proof**, in the
order the upcoming core paper introduces them. Each top-level
`stepN_*/` folder is a *container* — adding a new sub-experiment under
the same step never forces renames elsewhere.

For the precise theorem and proof strategy, read the root README §1–2
first.

## The folders, in proof order

| Folder | Proof step | Question it answers |
|---|---|---|
| [`step0_sanity/`](step0_sanity/) | §2. Setup sanity | Does $A_n^{\text{tree}}$ really vanish on each $\mathcal Z_r$? (the forward direction of Rodina's theorem) |
| [`step1_layer0_kill/`](step1_layer0_kill/) | §4 + §5. Step 1 (kill) | Which size-$(n{-}3)$ chord multisets does Layer-0 kill, and which escape? Plus a "dual" experiment isolating one zone's contribution. |
| [`step2_equate/`](step2_equate/) | §4 + §5. Step 2 (equate) | What equivalence classes does the bare↔special swap induce on Step-1 survivors? (The unitarity engine.) |
| [`step3_laurent/`](step3_laurent/) | §6. Step 3 (Laurent, n=7) | At $n = 7$, do all seven Step-1 survivors die at the next Laurent order? (Yes, all seven "fish".) |
| [`step4_laurent_block_analysis/`](step4_laurent_block_analysis/) | §6 + §7. Block-rule Laurent kill (n=8, n=9) | Does the depth-1 cascade rule extend? (n=8: 100/100 singletons. n=9: 94/113 singletons + 4 failures, but the 90-orbit cluster has rank 90 / nullity 0 — block-rule kills jointly.) |
| [`full_nullspace_verification/`](full_nullspace_verification/) | §5. End-to-end check | Does the full constraint system have nullspace dim $= 1$, equal to $A_n^{\text{tree}}$, at $n = 5, 6, 7$? |
| [`survivor_gallery/`](survivor_gallery/) | §5. Figures | Hand-curated PDFs of the $n = 8, 9$ non-triangulation survivor diagrams. |
| [`old/`](old/) | — | Superseded experiments tied to the older 15-dim / BCFW-induction approach. Kept for provenance. |

## Headline results (quick reference)

### Step-1 kill enumeration (`step1_layer0_kill/kill_enumeration/`)

| $n$ | total multisets | tri (Catalan) | non-tri survivors | non-tri kill % |
|---|---:|---:|---:|---:|
| 5 | 15 | 5 | **0** | 100% |
| 6 | 165 | 14 | **0** | 100% |
| 7 | 2 380 | 42 | **7** | 99.700% |
| 8 | 42 504 | 132 | **100** | 99.764% |
| 9 | 906 192 | 429 | **1 011** | 99.888% |

Triangulation counts match the Catalan numbers $C_{n-2}$ exactly. The
non-tri kill rate appears to grow with $n$ — the central conjecture is
that it goes to $1$ fast.

### Dual experiment (`step1_layer0_kill/dual_X13_never_special/`)

Excluding any single zone (default $\mathcal Z_{1,3}$), extra non-tri
survivors appear: 21 at $n=7$, 172 at $n=8$, 1 388 at $n=9$. **All
extras are non-triangulations** — no tree-amplitude diagram is ever
uniquely killed by a single zone, supporting the cyclic-symmetry
expectation.

### Step-2 equivalence classes (`step2_equate/equivalence_classes/`)

Conservative classes (only direct survivor↔survivor swaps) at $n=7$:
3 cyclic orbits, sizes $7\times\{1, 2, 4\}$. The remaining merging into
a single unitarity class is captured by the full-ansatz nullspace check.

### Laurent cascade (`step3_laurent/cascade_n7/`)

All 7 fish at $n=7$ die at order $S^{-3}$ (depth 1) in the Laurent
expansion, each via a 6-term fingerprint equation whose 5 cousins are
killed by Layer-0 at one neighbouring zone. Verified symbolically by
`cascade_kill_n7.py`; full trace in `results_cascade_n7.txt`.

### Full nullspace (`full_nullspace_verification/`)

At every $n \in \{5, 6, 7\}$, the full constraint system has nullspace
dimension $1$, spanned by $A_n^{\text{tree}}$. SVD gap $> 13$ orders of
magnitude.

## Quick run guide

```bash
python3 step0_sanity/verify_tree_zeros/verify_tree_zeros.py
python3 step1_layer0_kill/kill_enumeration/step1_uncaught.py 7
python3 step1_layer0_kill/dual_X13_never_special/step1_dual.py 7
python3 step2_equate/equivalence_classes/step2_classes.py 7
python3 step3_laurent/cascade_n7/cascade_kill_n7.py
python3 full_nullspace_verification/full_ansatz_n7.py
```

All scripts that produce figures write into a sibling `outputs/` folder
and (on macOS) auto-open the result.

## Soft cap on n

- $n \le 9$: easy (under a minute, $< 1$ MB PDFs).
- $n = 10$: slow ($\approx 10$ minutes, several thousand survivor diagrams).
- $n \ge 11$: hours; not recommended without optimization.
