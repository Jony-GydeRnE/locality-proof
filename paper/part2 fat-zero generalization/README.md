# Part II — the fat-zero (k-zero) generalisation

`fat_zero_part2.tex` → `fat_zero_part2.pdf` (6 pp).

**Core point.** Promote the hidden 1-zero to a **k-zero** (impose the hidden-zero relations on
*k consecutive boundaries*). The constraint surface $\mathcal Z^{(k)}$ has near block
`N={1..k+1}`, far block `F={k+2..n}`, length-k **enclosing chord** `X_{1,k+1}`, special
`X_{1,k+2}`, bare `X_{k+1,n}`. Rate-vector analysis (Theorem 1) gives the **complete**
classification of admissible asymptotic limits — every $\mathcal F$ comes from a rate vector
$(\alpha_\ell;\beta_i)\in\mathbb R^{(n-k-3)+(k-1)}$ via an explicit per-chord "finite iff" table.
The three named regimes:

- **(A)** all $\beta_i=0, \alpha_\ell=0$ — the 1-side;
- **(B)** all $\beta_i=1, \alpha_\ell=1$ — the (k+1)-side;
- **(C)** any other consistent rate vector (per-column slot picks + offset pairings).

Geometric shadow (Proposition 1): in the non-crossing case, survival forces a $(k{+}2)$-gon
completion of $N$, with **three** flavors — via the special, via the bare, or via an internal
**column pair** $\{X_{1m},X_{k+1,m}\}$ for any $m\in\{k{+}3,\dots,n{-}1\}$ (the column pair forces
$\alpha_m=0$ and $\alpha_m=1$ simultaneously).

**Build:** `python3 figures/make_figures.py` then `tectonic fat_zero_part2.tex`. Companion to
Part I (`final-paper-draft.pdf`, root) and Part III (`paper/part3 step2 K-zero relations/`).

**Source notes:** 20 (fat zeros), 21 (conjecture), 22 (P1/P2 dichotomy), 24-25 (A/B/C
classification), 26-27 (Step 2 — see Part III).

**Revision log.**
- v1 (2026-05-25): initial draft; Theorem 1 overclaimed exhaustiveness of A/B/C via a per-column
  three-slot rule that missed mixed-$\beta$ pairings; Prop 1 only invoked the crude limit.
- v2 (2026-05-30): rewrote Theorem 1 as the rate-vector classification (Eq. 7), upgraded Prop 1
  to all three $(k{+}2)$-gon completions including the column-pair (proof via rate-equation
  inconsistency). Survivors section now shows three families: special/bare, column pair, offset
  pair.
