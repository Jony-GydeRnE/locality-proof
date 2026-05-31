# Part II — the fat-zero (k-zero) generalisation

`fat_zero_part2.tex` → `fat_zero_part2.pdf` (5 pp).

**Core point — purely geometric.** Promote the hidden 1-zero to a **k-zero** (impose
hidden-zero on *k consecutive boundaries*). The constraint surface $\mathcal Z^{(k)}$
has near block `N={1..k+1}`, far block `F={k+2..n}`, length-$k$ **enclosing chord**
`X_{1,k+1}`, special `X_{1,k+2}`, bare `X_{k+1,n}`.

**Survival theorem (Theorem 1):** for every $k$ consecutive edges, a multiset $M$
survives iff at least one of:

- **(P1)** $M$ contains a chord crossing the enclosing chord, or
- **(P2)** $M$ extends the $(k{+}1)$-gon to a $(k{+}2)$-gon, in one of three flavours:
  - via the **special** $X_{1,k+2}$ (closes on the small side),
  - via the **bare** $X_{k+1,n}$ (closes on the large side),
  - via a **column pair** $\{X_{1,m}, X_{k+1,m}\}$ for some $m\in\{k{+}3,\dots,n{-}1\}$.

**Stacking.** The $(k{+}2)$-gon that closes a survival at the $k$-zero IS the
$(k{+}1)$-gon at the next $(k{+}1)$-zero. So the family of $k$-zero survival
conditions builds recursively: 2-gon → 3-gon → 4-gon → 5-gon → ...

**Multi-region P1.** When P1 must do the work, a single crossing isn't enough if
the diagram has multiple open regions — survival requires multiple chords crossing
the enclosing structure simultaneously, **one cherry per region** (note 22 page 6).

**Companion:** Part I (`paper/part1 1-zero asymptotic analysis/`) and Part III
(`paper/part3 step2 K-zero relations/`). Together with the end-to-end rank
experiment at `computations/full_constraint_system/`, this is "Paper 2" for
Rodina.

**Source notes:** 20 (fat zeros), 21 (the geometric survival statement; three
escapes), 22 (P1/P2 dichotomy + multi-region/cherry rule), 23 (case analysis of
two-step-removed survivors). Older draft revisions (with the rate-vector / α/β
framing) are superseded by this purely geometric version.

**Build:** `python3 figures/make_figures.py` then `tectonic fat_zero_part2.tex`.

**Sympy stress tests** at `tests/stress_tests.py` independently verify every
load-bearing claim (telescoping lemma, geometric survival = special/bare/column-pair,
binomial trace, per-column and offset Step-2 relations, triangulation escape).
All pass. The tests still use a rate-vector parametrisation internally as a
convenient scaffold for the symbolic checks; the paper itself does not.

**Revision log.**
- v1 (2026-05-25): first draft with per-column three-slot rule.
- v2 (2026-05-30): rewrote with a rate-vector classification (since dropped --- that
  notation was Giuseppe's bookkeeping and never made it into the notes).
- v3 (2026-05-31): rewrote around the geometric P1/P2 theorem + stacking +
  multi-region rule, matching notes 20-22 directly.
