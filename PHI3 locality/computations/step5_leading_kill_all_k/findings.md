# Findings — leading fat-zero kill at every $(k, i)$

**Verdict: World A confirmed at $n \le 9$, and the result is sharper than
the original hypothesis.** Locality at $n \le 9$ from cyclic hidden zeros
follows from the **leading fat-zero kill at $k=1$ and $k=2$ alone** — no
higher $k$, no Laurent cascade, no block rule, no Step-2 / Part III
relations needed.

## Headline numbers

| $n$ | total mset | triang | non-tri | step-1 survivors | killed by $k=2$ | survive ALL $(k,i)$ |
|----:|-----------:|-------:|--------:|-----------------:|----------------:|--------------------:|
| 5   | 15         | 5      | 10      | 0                | —               | **0** |
| 7   | 2 380      | 42     | 2 338   | 7                | **7**           | **0** |
| 8   | 42 504     | 132    | 42 372  | 100              | **100**         | **0** |
| 9   | 906 192    | 429    | 905 763 | 1 011            | **1 011**       | **0** |

**Every step-1 survivor's minimal-$k$ killing zone is $k=2$.** (Histogram
below.)

## The actual hypothesis was right; the conclusion is stronger

The user's gameplan asked: does *any* non-triangulation survive every
$(k,i)$ leading kill? Two possibilities:
- **World A**: none does.
- **World B**: some do.

Answer at $n\le 9$: **World A, and only $k=1$ and $k=2$ are needed.**

- $k=1$ kills 99.x% of non-triangulations directly (the Step-1 / Layer-0
  kill of Part I).
- $k=2$ kills *every single remaining survivor* — 100% of the 1 118
  step-1 survivors across $n=7,8,9$.
- $k\ge 3$ zones do additional work (most surviving multisets are
  killable at multiple $(k,i)$), but they are **never necessary** for
  the multisets that the $k=1$ analysis missed.

## Why this matters

Part I's all-of-$n=9$ argument required the **90-orbit block-rule cluster
matrix** (a $506\times 90$ rank-90 calculation) to handle the survivors
of the $k=1$ kill that didn't fall to a single-orbit Laurent cascade.
The block-rule machinery is **not needed if you apply the $k=2$ fat-zero
kill**. Every one of those 1 011 survivors dies under a single $k=2$
zone:

```
n=7  step-1 survivors (7)    : all killed by k=2  (none needs k>=3)
n=8  step-1 survivors (100)  : all killed by k=2  (none needs k>=3)
n=9  step-1 survivors (1011) : all killed by k=2  (none needs k>=3)
```

So the *leading* fat-zero kill — Part II's geometric P1/P2 criterion at
$(k, i)$ — applied at $k\in\{1, 2\}$ across all cyclic rotations is
sufficient for locality at $n \le 9$. Part III (Step-2 coefficient
equations) is *only* needed for unitarity (equating triangulation
coefficients), not for locality.

## Histogram of minimal killing $(k, i)$ at $n=9$

For each of the 1 011 step-1 non-triangulation survivors, the
lexicographically minimal $(k, i)$ that kills it (k=2 for every one):

| $i$ | count |
|----:|------:|
| 1   | 337   |
| 2   | 296   |
| 3   | 108   |
| 4   | 105   |
| 5   | 89    |
| 6   | 38    |
| 7   | 26    |
| 8   | 8     |
| 9   | 4     |

Most fall to a $k=2$ zone based at $i=1$ or $i=2$. The skew is partly an
artifact of the lexicographic minimisation (early-$i$ wins ties); each
survivor is typically killable by **many** $(k,i)$ zones simultaneously.

## What was used

**PRIMARY check:** rate-equation feasibility (sympy `linsolve`). A
multiset is killable at $(k, i)$ iff the linear system on the rate
variables $(\alpha_\ell, \beta_p)$ imposed by "each chord finite" has a
solution. This is the rigorous criterion derived in Part II.

**Notes on the geometric criterion.** Part II's v4 Theorem 1 states
"survive iff P1 (a chord crosses the enclosing chord) or P2 ((k+2)-gon
completion)". The "P1 ⇒ survives" half is too weak: a single crossing
chord can typically be held finite by a careful rate choice. The
*correct* survival conditions are the chord-pair incompatibilities of
v4 §5 (special / bare / column pair / offset pair), **plus higher
multi-chord conflicts** that arise from crossings paired with offsets.
The rate-equation feasibility check captures all of them; the geometric
V1–V4 alone is sufficient but not necessary.

This script uses the rate-equation check throughout, so the World A
verdict is unaffected.

**Self-check passed.** Restricted to $k=1$, the code reproduces the
known step-1 survivor counts: 7 / 100 / 1 011 at $n=7$ / 8 / 9. Match
to the existing Part I enumeration (`computations/step1_layer0_kill/`).

## Outputs

```
outputs/
  n5_findings.md, n5_phase0.json, n5_phase1.json
  n7_findings.md, n7_phase0.json, n7_phase1.json
  n8_findings.md, n8_phase0.json, n8_phase1.json
  n9_findings.md, n9_phase0.json, n9_phase1.json
  summary.json
```

The `n*_phase0.json` files contain, for each step-1 non-triangulation
survivor, the full list of $(k, i)$ zones that kill it plus the
lexicographically minimal one. This is the raw material for the
structural "band-cover lemma" that would close the all-$n$ proof.

## What's next (the structural lemma)

The empirical observation "every step-1 survivor at $n \le 9$ dies at
$k=2$" begs to be promoted to a theorem. Concretely:

> **Band-cover conjecture.** Every non-triangulation multiset $M$ of
> the $n$-gon contains either (a) the special / bare / column-pair /
> offset-pair of *some* $k$-zone, $k \in \{1, 2\}$, or (b) a chord that
> participates in a rate-equation incompatibility with another chord
> of $M$ at some $(k, i)$, $k \in \{1, 2\}$.

If true, locality at every $n$ follows from the $k=1$ and $k=2$ leading
kills alone — no induction over $n$, no residue argument, no Step-2.

The next step is to derive this directly from the geometry of the
$n$-gon, using the chord-crossing structure of non-triangulations. This
would supersede Part I's per-$n$ Laurent cascade / block-rule and
produce a much cleaner uniform-in-$n$ statement of locality.

## Reproducing

```bash
cd computations/step5_leading_kill_all_k
python3 leading_kill_all_k.py --ns 7 8                    # ~12 s
python3 leading_kill_all_k.py --ns 9 --skip-phase1 \      # ~2 min (phase 0 only)
                              --skip-selfcheck            #
# For the full Phase-1 enumeration at n=9 (~110 s), drop --skip-phase1.
```

Outputs land in `outputs/`. Deterministic; no randomness.
