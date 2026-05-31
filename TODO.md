# TODO — pickup notes for future sessions

Last touched: 2026-05-31.

## Current state of the repo

- **Three-part write-up is live.** All in `paper/`:
  - **Part I** — `paper/part1 1-zero asymptotic analysis/one_zero_part1.pdf` (7 pp).
    Locality from the cyclic 1-zero, proven through $n=9$ (Step-1 kill +
    triangulation escape + asymptotic completeness $S(n)/T(n)\le 8/n$ +
    $n=7$ fish + $n=9$ 90-orbit block rule).
  - **Part II** — `paper/part2 fat-zero generalization/fat_zero_part2.pdf` (6 pp).
    The $k$-zero generalisation. Rate-vector classification of admissible
    asymptotic limits → criteria $A$, $B$, $C$. Geometric shadow: three
    $(k{+}2)$-gon completions of the near block (special, bare, column pair).
  - **Part III** — `paper/part3 step2 K-zero relations/step2_part3.pdf` (6 pp).
    Step-2 coefficient relations on each $k$-zero. Universal binomial trace
    $\sum_j (-1)^j a_{s^{m-j}b^j}=0$ + per-column + per-offset relations.
    Worked for $k=1,2,3$.
- **Symbolic stress tests** at `paper/part2 fat-zero generalization/tests/stress_tests.py`
  pass on every load-bearing claim of Parts II/III: telescoping lemma,
  chord-rate table, Prop 1 across 17,602 multisets, binomial trace,
  per-column / offset relations, triangulation escape across 400
  Catalan-many triangulations.
- **All older paper drafts** (step1 / step2 Laurent / step3 block rule /
  Details missing / elementary-geometric) live in
  `old stuff/paper-older-works/` (gitignored).

## The open question

**All-$n$ locality** from the union of $k$-zero kills (Part II) plus
Step-2 relations (Part III). Concretely: does the combined linear system
on the multiset coefficients $a_M$ have cokernel exactly on the
triangulation coefficients, uniformly in $n$? Empirical data through
$n=9$ is consistent.

## Concrete next steps

- Promote the Theorem 2 (per-$m$ structure) sketch in Part III §6 to a
  full induction; this is where the "polynomial growth" $O(n^3)$ count
  is currently back-of-envelope.
- Tighten the mixed regime $C$ statement in Part II §4 — the per-column
  three-slot picture (Fig 4) is canonical, but the full $C$ family from
  the rate vector includes additional pairings (verified in T2 of
  stress_tests). Spell this out in a corollary.
- Extend the stress-test suite to $n=10$, $k=4$ (currently $n\le 8$,
  $k\le 4$). Cost: T2 enumeration grows quickly; may need filtered
  enumeration.
- Optional: clean source-of-Part-I — `one_zero_part1.pdf` is currently
  the built artifact only. If the `.tex` source is wanted in
  `paper/part1 .../`, recover from `old stuff/paper-older-works/step1
  Kill technique and statistics/locality_unitarity_v6_DRAFT.tex`.

## Computations housekeeping

- The `computations/` folder is current; all subfolders back Part I.
  `step1_layer0_kill/short_edge_property/` was untracked before this
  reorg — included in the next commit.
- $n=10$ block-rule extension is still untried.

## What NOT to do

- Don't touch `notes/`, `existing literature/`, `old stuff/` contents.
- Don't rebuild the Step-2 flip graph or the $n=9$ cluster matrix —
  those are done and recorded; pull from `outputs/`.
- Don't try to "kill" the 4 fan-class triangulation components at
  $n=9$ — those are LOCAL and survive by design.
