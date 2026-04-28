# §5 — Step-1 kill statistics and the asymptotic conjecture (Item 4)

**Backs paper §5 (numerical evidence) and Item 4 of the master plan.**

This satellite documents the *global* kill rate of the Step-1
mechanism: at each $n$, what fraction of all non-triangulation
size-$(n-3)$ multisets does Step-1 catch, and how does that fraction
grow with $n$?

It is the global counterpart to
[`../step-one-doesnt-kill-triangulations/`](../step-one-doesnt-kill-triangulations/),
which proves the *local* side (no triangulation is ever Step-1 killable).

## Headline numbers (from `computations/step1_layer0_kill/kill_enumeration/`)

| $n$ | non-tri total | non-tri survivors | non-tri kill % |
|---|---:|---:|---:|
| 5 | 10 | 0 | 100% |
| 6 | 151 | 0 | 100% |
| 7 | 2 338 | 7 | 99.700% |
| 8 | 42 372 | 100 | 99.764% |
| 9 | 905 763 | 1 011 | 99.888% |

Kill rate is monotonically increasing in $n$. Survivor fraction
$S(n)/T(n)$ roughly halves per increment of $n$.

## Files

- **`step1_statistics.tex`** — the satellite paper.
- **`step1_statistics.pdf`** — compiled (recompile after editing).

## What it argues

1. The Step-1 mechanism handles the overwhelming majority of
   non-triangulation coefficients, with kill rate $\to 1$.
2. The asymptotic conjecture is stated as Conjecture 1 of the paper:
   $S(n)/T(n) \le C \cdot 2^{-n}$ (exponential decay favoured by data,
   polynomial decay not ruled out by $n \le 9$).
3. A proof would require characterising the cyclic-orbit family of
   survivors and bounding its growth.
4. **Important:** this conjecture is *not* load-bearing for the
   all-$n$ locality theorem. What is load-bearing is the per-survivor
   Laurent cascade (Item 10), verified at $n = 7$ and being verified
   at $n = 8$ in `computations/step3_laurent/cascade_n8/`.

## Status

- $n = 7$: 7 survivors, all killed by depth-1 Laurent cascade
  (`paper/Laurent series for hard kills/`).
- $n = 8$: **100 / 100 survivors killed by depth-1 Laurent cascade**
  (verified symbolically in 76 min, see
  `computations/step3_laurent/cascade_n8/`). 90 survivors die at
  Laurent order $X_*^{-3}$, 10 at order $X_*^{-4}$. Cumulative
  Step-1 + Step-3 kill rate on $n=8$ non-triangulations: **100%**.
- $n = 9$: 1011 survivors; cascade verification not yet attempted
  (will require adapting the $n=8$ driver and longer compute time).
- Conjecture 1 (asymptotic decay of $S/T$): **central open** — code
  provides data and motivates the bound but cannot prove it.
