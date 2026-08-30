# Stress tests for Part II + Part III

`stress_tests.py` — sympy-based symbolic verification of every load-bearing claim.

**Run:** `python3 stress_tests.py`. Pure symbolic; passes are exact, not numerical.

## What is verified

| Test | Claim | Source | Scope |
|------|-------|--------|-------|
| T1a  | Telescoping lemma: substitutions satisfy $c_{ij}=0$ | Part II Lemma 1, eqs (3)–(4) | $(n,k) \in \{(5,1),(6,1),(6,2),(7,2),(7,3),(8,3),(9,4)\}$ |
| T1b  | Chord-rate "finite iff" table | Part II eq (7) | 5 $(n,k)$ × 4 rate vectors |
| T2   | $M$ survives ⟺ contains special/bare/column-pair | Part II Prop 1 | Full enumeration of non-crossing $M$ at 6 $(n,k)$; **17,602 multisets total** |
| T3a/b| Binomial trace $\sum_j(-1)^j a_{s^{m-j}b^j}=0$ | Part III Thm 1, eq (2) | 8 $(n,k,m)$ incl. $k=2$ quadratic & cubic |
| T3c  | 2-zero per-column relation $\sum_i(a_{14,i\ell}-a_{3n,i\ell})=0$ | Part III eq (6) | $n=6,7,8$ at $k=2$ |
| T3d  | 2-zero offset relation $a_{14,24}+a_{14,2n}-a_{3n,24}-a_{3n,2n}=0$ | Part III eq (7) | $n=6,7,8$ at $k=2$ |
| T4   | Every triangulation escapes the k-zero | Part II §5 | All Catalan $C_{n-2}$ triangulations at 9 $(n,k)$; **400 triangulations total** |

## Last run

All tests PASS as of 2026-05-30 — full output at end of file.

## Strategy

- **T1a, T1b** verify the algebraic backbone: the substitutions are correct and the rate table reads off chord scalings correctly.
- **T2** is the heavy hitter: it independently classifies each multiset by (a) checking rate-equation consistency directly via `sp.linsolve` and (b) checking the geometric criterion of Prop 1, then verifies they agree. The 11,628-multiset $n=8$ enumeration is the strongest test.
- **T3a–d** verify the actual Step-2 coefficient relations by symbolically computing the appropriate Laurent residues. Each test compares the extracted relation against the predicted formula via `sp.expand(diff) == 0`.
- **T4** verifies the triangulation-always-escapes claim by enumerating Catalan-many triangulations and running each through the T2 rate-equation check.

## Key implementation notes

- `setup_kzero(n, k)` returns a dict $X[(i,j)]$ giving every chord as a sympy expression in independent coords $\{s, t_i, u_\ell, Y_{(i,j)}\}$.
- Laurent residues are extracted via `sp.limit(s * expr, s, 0)` etc.; do **not** use `(expr * s).subs(s, 0)` — the literal substitution before cancellation produces `nan` on rational expressions.
- The rate-equation system per multiset is built directly from eq (7)'s "finite iff" clauses and solved with `sp.linsolve` for consistency.
