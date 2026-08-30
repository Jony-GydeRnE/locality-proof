# `full_constraint_system/` — the union of all k-zero constraints

**One-line takeaway.** At $n = 5, 6, 7$, the rank of the full constraint
matrix (vanishing of $B_n$ on every $k$-zero, every cyclic rotation)
equals $T(n) - 1$, and the unique nullspace direction is exactly the
triangulation indicator. So the full $k$-zero system **closes** at
$n = 7$: it forces every non-triangulation coefficient to zero and
equates every triangulation coefficient to a single common scalar
— locality *and* unitarity simultaneously, from the union of the
$k$-zero constraints alone, **strictly stronger than Part I's $k=1$
result.**

## What this answers

Parts II and III each end with the same open question:

> Does the union of the $k$-zero kills (Part II) plus the Step-2
> coefficient relations (Part III) close to a 1-dimensional cokernel
> spanned by the triangulation indicator?

This folder turns that into a sparse linear-algebra computation and
runs it at small $n$.

## The unification

Both Parts II and III are projections of one bigger object:

> $B_n |_{\mathcal Z^{(k)}_r} \;=\; 0$ as a Laurent polynomial in the free
> coordinates of $\mathcal Z^{(k)}_r$.

Each monomial coefficient of that vanishing is a linear functional on
the multiset coefficients $a_M$ that must equal zero. Part II picks off
some of these (the rate-vector kills); Part III picks off others (the
$X\!\to\!0$ Step-2 relations). Neither is exhaustive on its own.

This script computes the rank of the *combined* system directly, by
sampling random rational points on each $\mathcal Z^{(k)}_r$ — each
sampled point gives one linear functional
$\sum_M (1/\prod_{c\in M}X_c)\,a_M = 0$, and the rank of the matrix of
these functionals over $\mathbb F_p$ tells you the dimension of the
*total* constraint space at this $n$.

## Method

For each $(k, r)$ with $k \in \{1,\dots,n{-}3\}$ and $r \in \{1,\dots,n\}$:

1. Parametrise $\mathcal Z^{(k)}_r$ by the free coords (special,
   offsets, column-1, born-free). All canonical chord labels are
   cyclically shifted by $r$.
2. Sample random rational values for the free coords.
3. Compute every chord $X_{i,j}$ via the cyclic-shifted Part II Lemma 1
   substitutions.
4. Build the constraint row $\text{row}[\text{idx}_M] = (\prod_{c\in M} X_c)^{-1} \bmod p$.
5. Add to an incremental row-reduced echelon over $\mathbb F_p$
   ($p = 998244353$).

Sampling continues until the rank either reaches $T(n) - 1$ (the
expected value) or stops growing across multiple sweeps over the
$(k, r)$ pairs.

The conjecture is **rank $= T(n) - 1$ with nullspace $=$ span of the
triangulation indicator** (which is verified directly by checking
that the triangulation indicator vector lies in the row space of every
pivot row, mod $p$).

## Results

| $n$ | $T(n)$ | $C_{n-2}$ | $(k,r)$ pairs | samples | rank | nullspace dim | triang indicator $=$ kernel | time |
|---:|---:|---:|---:|---:|---:|---:|:---:|---:|
| 5 | 15    | 5  | 10 | 14    | **14** | **1** | ✓ | 0.0 s |
| 6 | 165   | 14 | 18 | 164   | **164** | **1** | ✓ | 0.2 s |
| 7 | 2 380 | 42 | 28 | 2 379 | **2 379** | **1** | ✓ | 70 s |

At every $n$ tested, the **unique kernel direction**
$\ker M$ has exactly $C_{n-2}$ nonzero entries, all on triangulations,
all equal to the same value mod $p$. That is, **the only solution to
the full vanishing constraints is "$a_M = c$ for every triangulation
$M$, $a_M = 0$ otherwise"** — locality and unitarity in one stroke.

Per-$n$ markdown reports: `outputs/n{5,6,7}_results.md`. JSON
summaries: `outputs/n{5,6,7}_results.json` and `outputs/summary.json`.

## Significance

- **Strictly stronger than Part I.** Part I proves locality at $n \le 9$
  from the $k = 1$ zeros (with $n = 9$ needing the 90-orbit block-rule
  cluster). Here we prove locality + unitarity *together*, at every
  $n \le 7$, from the union of $k = 1, \dots, n-3$ zeros — no Laurent
  cascade, no block rule needed. The fat zeros do all the work.
- **Empirical confirmation that Parts II/III together close the
  system.** Both papers ended saying "the open question is whether the
  combined system closes." Now we know: at every $n$ tested, it does.
- **A clean target for the all-$n$ proof.** With the empirical answer
  in hand, the structural theorem becomes "for every $n$, the
  $T(n) \times T(n)$ matrix of all (chord-evaluation, $(k,r)$) pairs
  has rank $T(n) - 1$ with the triangulation indicator in its
  kernel." This is a much sharper open question than "does Part I's
  block-rule cascade keep working."

## $n = 8$ and beyond

$n = 8$ has $T(8) = 42{,}504$ unknowns. With the current
$O(\text{rank}^2 \cdot T)$ row-reduction this would take ~hours-to-days
in Python; it was attempted and killed. Next steps to push to $n = 8, 9$
are listed in the parent `TODO.md`:

- Use FLINT / sage / a C extension for fast modular linear algebra.
- Switch to a randomized rank algorithm (Wiedemann / block Lanczos).
- Or block-decompose by cyclic-orbit and reduce within orbits first.

## Files

- `build_full_system.py` — the script (commented).
- `outputs/n{5,6,7}_results.md` — human-readable per-$n$ verdicts.
- `outputs/n{5,6,7}_results.json` — machine-readable summaries, including
  the nullspace-basis structure (nonzero counts split into triangulation
  vs non-triangulation, and the "matches_triang_indicator" flag).
- `outputs/summary.json` — combined summary across all $n$.

## Reproducing

```bash
cd computations/full_constraint_system
python3 build_full_system.py 5 6 7
```

Deterministic seed (42 by default). $n=7$ takes about a minute on a
typical laptop; $n=5,6$ are instant. Outputs land in `outputs/`.

## Caveat: probabilistic rank

The rank is computed by sampling rational evaluations and reducing mod
a large prime. This gives the rank exactly (with overwhelming
probability) as long as no chord value happens to be zero mod $p$
(handled by regenerating points) and as long as the prime doesn't
divide any non-zero minor of the constraint matrix (extremely unlikely
for $p \approx 10^9$). For absolute certainty one would re-run with a
second prime and confirm agreement; in this folder we use one prime
($p = 998244353$).
