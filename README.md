# Locality and Unitarity of Tr(φ³) Tree Amplitudes from Hidden Zeros
By Jony V

> A working repository building toward an algebraic proof that the cyclic
> hidden 1-zero conditions uniquely determine the Tr(φ³) tree amplitudes —
> locality and unitarity emerge as *consequences*, not assumptions.
>
> The repo doubles as (i) an entry point into Rodina's paper
> ([arXiv:2406.04234](https://arxiv.org/abs/2406.04234)) — turning the
> implicit steps of his proof into self-contained arguments aimed at
> upper-level undergraduates — and (ii) an extension of his hidden-zero
> theorem into a direct proof of locality verified through $n = 9$.

This README is the entry point. Read it top to bottom: it states the problem
precisely, gives the proof strategy in plain words, summarizes the current
status, and walks through every folder in the repo and how it supports the
conjecture.

---

## 📌 Reader's guide — where each result lives

Every claim in the accompanying letter is backed by a verifiable file in
this repo. Use this table to jump directly from a claim to the file that
establishes it.

| Result | File / folder |
|---|---|
| **Part I — locality from the 1-zero, through $n=9$.** Step-1 asymptotic kill (Theorem 1), triangulation escape (Theorem 2), asymptotic completeness $S(n)/T(n)\le 8/n+O(n^{-2})$, the $n=7$ fish (Theorem 3), the $n=9$ 90-orbit block rule. | [`paper/part1 1-zero asymptotic analysis/one_zero_part1.pdf`](paper/part1%201-zero%20asymptotic%20analysis/) |
| **Part II — the fat-zero generalization.** Promote the 1-zero to a $k$-zero. Rate-vector classification of admissible asymptotic limits gives criteria $A$, $B$, $C$; the non-crossing survival builds a $(k{+}2)$-gon over the length-$k$ enclosing chord (special / bare / column-pair). | [`paper/part2 fat-zero generalization/fat_zero_part2.pdf`](paper/part2%20fat-zero%20generalization/) |
| **Part III — Step-2 $k$-zero relations.** Universal binomial trace $\sum_j(-1)^j a_{s^{m-j}b^j}=0$ at every $k,n,m$; per-column and per-offset relations; the genuinely new $k\ge 2$ identities like $a_{14,14}+a_{3n,3n}-a_{14,3n}=0$. | [`paper/part3 step2 K-zero relations/step2_part3.pdf`](paper/part3%20step2%20K-zero%20relations/) |
| **Symbolic stress tests.** Sympy verification of every load-bearing claim in Parts II/III: telescoping lemma, chord-rate table, Prop 1 over 17,602 multisets, binomial trace, per-column / offset relations, triangulation escape over 400 Catalan-many triangulations. | [`paper/part2 fat-zero generalization/tests/stress_tests.py`](paper/part2%20fat-zero%20generalization/tests/) |
| **Step-1 kill statistics across $n$.** 100% at $n \le 6$, 99.7–99.9% at $n = 7, 8, 9$. | Part I §4 + `computations/step1_layer0_kill/kill_enumeration/` |
| **$n=9$ block-rule per-orbit audit.** All 113 orbit reps accounted for: 23 single-orbit cascades + 90-orbit block rule (rank 90 / nullity 0). | `computations/step4_laurent_block_analysis/n9/outputs/n9_locality_status.md` |
| **Full nullspace dim $= 1$ at $n = 5, 6, 7$.** End-to-end SVD verification of the complete constraint system. | `computations/full_nullspace_verification/` |
| **Reproducibility.** Every script that produced a number in this repo is committed alongside its raw output and a folder README. | `computations/` (top-level guide); each subfolder has its own README |
| **Working notes (~70 PDFs).** Notes 12, 16–27 are the load-bearing notes for Parts II/III; notes 13–17 for Part I. | `notes/` |

The status table below (§3) gives the same picture in tabular form, and
the guided tour (§4) walks through how to read the repo end-to-end.

---

## 1. The problem we want to prove

Label the legs of an n-point amplitude cyclically by 1,…,n. Write the
planar Mandelstam invariants $X_{ij} = (p_i+\dots+p_{j-1})^2$ for the
$d = n(n-3)/2$ chords of the n-gon. Let

$$B \;=\; \sum_{w}\frac{a_w}{\prod_{c\in w} X_c}$$

be the most general rational function of mass dimension $-2(n-3)$ with
at most simple poles on planar invariants ($w$ runs over size-$(n-3)$
multisets of chord indices). Define the cyclic hidden 1-zero locus

$$\mathcal{Z}_r \;:=\; \{c_{r,r+2}=c_{r,r+3}=\dots=c_{r,r+n-2}=0\},$$

the codimension-$(n-3)$ subspace where the $n-3$ non-planar Mandelstams
incident to vertex $r$ vanish. Rodina (arXiv:2406.04234) proved that the
Tr(φ³) tree amplitude vanishes there: $A_n^{\text{tree}}\big|_{\mathcal Z_r} \equiv 0$
for every $r$.

**Theorem (target — what the upcoming core paper will prove).**
*If $B\big|_{\mathcal Z_r} \equiv 0$ for every $r \in \{1,\dots,n\}$, then*
$$B \;=\; c \cdot A_n^{\text{tree}}$$
*for some scalar $c$.*

Locality (only Feynman / triangulation poles survive) and unitarity (all
triangulation coefficients are equal) emerge **simultaneously** as
consequences of the $n$ cyclic 1-zero constraints — not as separate
inputs.

---

## 2. The proof strategy in one paragraph

At each zone $\mathcal Z_r$, the master substitution rewrites the $n-3$
chords incident to vertex $r{+}1$ as $\mathbb Q$-linear combinations of
the remaining (free) chords $F_r$ plus one *special* chord
$X_* := X_{r,r+2}$. Three layers of progressively finer arguments peel off
the coefficients of $B$:

- **Step 1 (Layer-0 kill).** Send $X_*$ and every chord outside $F_r$ to
  $\infty$, holding $F_r$ free. Surviving terms of $B|_{\mathcal Z_r}$ are
  Laurent monomials in $F_r$; distinct multisets in $F_r$ give distinct
  monomials, so each coefficient $a_M$ with $M \subseteq F_r$ is
  independently forced to $0$.
- **Step 2 (equate).** Keep $X_*$ in the surviving set. The bare
  relation $X_{r+1,r-1} = -X_*$ glues pairs of multisets that differ
  only by this swap, producing equalities $a_M = a_{M'}$. Cyclic shifts
  of these equalities form a flip-graph on the $C_{n-2}$ triangulations;
  flip-connectivity collapses all triangulation coefficients to one
  common value.
- **Step 3 (Laurent).** Some non-triangulation multisets at $n\ge 7$
  escape both Step 1 and Step 2 (the "fish" at $n=7$). Expanding
  $1/(X_{\text{companion}}-X_*)$ as a Laurent series in $1/X_*$, one
  cascades to higher Laurent orders to kill each escapee using
  already-established Step-1 kills at one neighbouring zone.

**Item 2 (formerly an important thing to verify) is now closed.**
Theorem 2 of Part I (triangulations escape) proves the strictly stronger
statement that *no triangulation $T$ is in $\mathcal K_{r,r+2}$ at any zone
$r$* — using only the polygon edge $(r, r{+}1)$ and the unique triangle of
$T$ that contains it.

The **central question** is the *all-$n$ uniformity* of the
Laurent cascade (Items 7 + 10 below):

> **Conjecture (Item 10).** For every $n \ge 7$, every Step-1 survivor
> dies via a depth-1 Laurent cascade of the same shape as the $n=7$ fish
> kills — *one Laurent step deeper at one zone, with all prerequisites
> of layer-0 type at one neighbouring zone*.

If true, every non-triangulation coefficient $a_M$ is forced to zero,
and combined with Step-2 flip-graph connectivity, the surviving space
is one-dimensional and spanned by $A_n^{\text{tree}}$. The depth-1
phenomenon is verified at $n=7$ symbolically; the $n=8$ verification
is in progress (`computations/step3_laurent/cascade_n8/`). At $n=9$
the survivor fraction is already 99.889%, and the rate appears to grow
fast with $n$ (Item 4).

---

## 3. Status by n

| $n$ | Step 1 alone | Step 1 + 2 | Step 1 + 2 + 3 | Full nullspace |
|---|---|---|---|---|
| 4 | trivial | ✓ | ✓ | dim 1 |
| 5 | all 10 non-locals killed | ✓ | ✓ | **dim 1 (verified)** |
| 6 | 129 / 151 killed | ✓ remaining 22 die by cyclic Step 2 | ✓ | **dim 1 (verified)** |
| 7 | 7 fish escape | 7 fish still escape | ✓ all 7 die by Laurent cascade | **dim 1 (verified)** |
| 8 | 100 escape | … | ✓ **all 100 die by depth-1 Laurent cascade** | — |
| 9 | 1 011 escape (out of 906 192 non-tri multisets) | not relevant for locality at this $n$ (unitarity handled separately by d-subset) | **all 1 011 survivors die: 23 single-orbit depth-1 cascades + 90-orbit block-rule cluster (rank 90, nullity 0)** | — |

> **Bottom line: locality is proven at $n \le 9$ from the cyclic 1-zeros
> alone** (Part I). Unitarity is recovered through the Step-2 D-subset
> mechanism (Part III, §3, and note 18).

At $n = 9$ Step 1 alone already kills **99.888%** of all non-triangulation
multisets. The hope is to prove the survivor fraction drops to $0$ fast as
$n$ grows, and that the residual "hard kills" can be handled uniformly by
the Laurent cascade.

---

## 4. How to read this repo (the guided tour)

Follow this order. Each section ends with a pointer to the relevant file
or subfolder. General notes are found in the Notes. This content takes one from trace phi 3 lagrangian and Feynman diagrams as triangulations to the the detailed proof that unitarity emerges at all n via hidden zeroes, following Rodina's paper https://arxiv.org/abs/2406.04234.  

### 4.1 Read Part I — locality from the 1-zero through $n=9$
→ [`paper/part1 1-zero asymptotic analysis/one_zero_part1.pdf`](paper/part1%201-zero%20asymptotic%20analysis/)

The Step-1 asymptotic kill, the triangulation-escape theorem, the
asymptotic completeness bound $S(n)/T(n)\le 8/n+O(n^{-2})$, the
worked $n=7$ "fish" via a depth-1 Laurent cascade, and the $n=9$
block-rule cluster argument (rank 90 / nullity 0 on the 90-orbit
cluster).

### 4.2 Read Part II — the fat-zero generalisation
→ [`paper/part2 fat-zero generalization/fat_zero_part2.pdf`](paper/part2%20fat-zero%20generalization/)

Promote the cyclic hidden 1-zero to a $k$-zero (impose hidden-zero on $k$
consecutive boundaries). Telescoping lemma + rate-vector classification
gives the criteria $A$, $B$, $C$. The non-crossing geometric shadow:
build a $(k{+}2)$-gon over the length-$k$ enclosing chord in three
flavours — via the special, the bare, or an internal column pair.

### 4.3 Read Part III — Step-2 $k$-zero coefficient relations
→ [`paper/part3 step2 K-zero relations/step2_part3.pdf`](paper/part3%20step2%20K-zero%20relations/)

Boundary-propagator grading + D-subsets. Universal binomial trace
$\sum_j(-1)^j a_{s^{m-j}b^j}=0$ at every $(k,n,m)$, plus per-column and
per-offset relations. At $k=1$, $m=1$ this recovers $a_{13}=a_{2n}$;
at $k\ge 2$, $m\ge 2$ it produces genuinely new relations such as
$a_{14,14}+a_{3n,3n}-a_{14,3n}=0$.

### 4.4 (Optional) Verify the claims symbolically
→ [`paper/part2 fat-zero generalization/tests/stress_tests.py`](paper/part2%20fat-zero%20generalization/tests/)

A sympy stress-test suite that symbolically checks every load-bearing
claim in Parts II/III: telescoping lemma, chord-rate table, Prop 1 by
exhaustive enumeration (17,602 multisets), binomial trace at $m=1,2,3$,
per-column and offset relations, triangulation escape across 400
Catalan-many triangulations.

### 4.5 Run the experiments
For numerical verification, survivor enumeration, and figures:

→ `computations/` — see `computations/README.md`

The folders are numbered to match the proof's logical flow:
`step0_sanity` → `step1_layer0_kill` → `step2_equate` → `step3_laurent`,
plus end-to-end checks (`full_nullspace_verification`) and figures
(`survivor_gallery`).

### 4.6 Browse the chronological notes
The complete working notebook (~70 PDFs):

→ `notes/`

Notes 1–27 are the topical PDFs (#1 = Feynman diagrams & triangulations,
#27 = Step-2 K-zero relations). Load-bearing for Parts II/III: notes 12,
16–27. Load-bearing for Part I: notes 13–17. The rest are dated session
notes from May 2025 onwards.

### 4.7 Existing literature
Background papers referenced throughout:

→ `existing literature/`

Includes Rodina (arXiv:2406.04234, the foundational paper this work
extends), Arkani-Hamed et al. (arXiv:2312.16282), Gonzales–Ward, and the
Feynman-diagram interpretation paper.

---

## 5. Repository layout

```
Rodina-locality-proof/
│
├── README.md                                ← you are here
│
├── paper/                                   ← three-part write-up
│   ├── README.md                              guided tour of paper/
│   ├── part1 1-zero asymptotic analysis/      Part I  — 1-zero kill, through n = 9
│   ├── part2 fat-zero generalization/         Part II — k-zero kill (A/B/C criteria) + tests/
│   └── part3 step2 K-zero relations/          Part III — Step-2 coefficient relations
│
├── computations/                            ← experiments backing each proof step
│   ├── README.md                              guided tour of computations/
│   ├── step0_sanity/                          §2 — A_n^tree vanishes on each Z_r
│   ├── step1_layer0_kill/                     §4 — Step-1 kill enumeration + dual
│   ├── step2_equate/                          §4 — Step-2 equivalence + flip-graph (n=9)
│   ├── step3_laurent/                         §6 — Laurent cascade at n=7
│   ├── step4_laurent_block_analysis/          §6 + §7 — n=8 singleton + n=9 block-rule
│   ├── full_nullspace_verification/           §7 — full-ansatz nullspace = 1
│   ├── survivor_gallery/                      §7 — n=8, n=9 survivor figures
│   └── old/                                   superseded experiments (15-dim approach)
│
├── notes/                                   ← ~50 chronological PDF notes
│
├── existing literature/                     ← Rodina, Arkani-Hamed, …
│
└── old stuff/                               ← historical drafts, parallel approaches
```

---

## 6. Open items (the master plan)

| # | Item | Status |
|---|---|---|
| 1 | Foundational lemma $B = 0 \Rightarrow B_i = 0$ | **closed** — note 12; underwrites the boundary-propagator grading of Part III |
| 2 | No triangulation is killed by Step 1 at any zone | **closed** — Part I Theorem 2; Part II Prop 1 generalises to every $k$-zero (verified for $n\le 8$ in `paper/part2/tests/stress_tests.py` T4) |
| 3 | Use Item 2 to identify exactly which indices Step 1 kills | now mechanical (Items 1 + 2 closed) |
| 4 | Survivor fraction $\to 0$ fast as $n$ grows | data through $n=9$ (Part I §4); asymptotic conjecture open |
| 5 | Numerical verification at $n = 7, 8, 9$ | done — `computations/step3_laurent/` and `computations/step4_laurent_block_analysis/` |
| 6 | Dual experiment ($X_{13}$ never special) | done — `computations/step1_layer0_kill/dual_X13_never_special/` |
| 7 | **Locality at $n \le 9$ from the 1-zero alone** | **closed** — Part I; Step-1 + 23 single-orbit cascades + block-rule 90-orbit cluster kill ⟹ every non-triangulation coefficient is zero. |
| 8 | **$k$-zero classification of kills** | **closed** — Part II Theorem 1 (rate-vector classification of admissible finite-sets); verified symbolically over 17,602 multisets at $n\le 8$ |
| 9 | **Step-2 $k$-zero relations** | **closed for $k=1,2,3$** — Part III; binomial trace + per-column + per-offset relations verified symbolically |
| 10 | All-$n$ locality from the union of $k$-zero kills + Step-2 relations | **open** — central conjecture; the question is whether the combined linear system has cokernel exactly on the triangulation coefficients |
| 11 | Triangulation $c$-values collapse to a single common scalar (UNITARITY) | recovered from Part III's $m=1$ binomial trace $a_{X_{1,k+2}}=a_{X_{k+1,n}}$ at every $k$ and rotation; the resulting flip-graph collapses triangulation classes |

The updated dependency chain is:

**Items 1 + 2 (closed) + Items 7–9 (closed at $n \le 9$ and $k \le 3$) ⟹ Item 10 (all-$n$ locality) + Item 11 (unitarity) ⟹ locality + unitarity theorem.**

---

## Acknowledgements

The hidden-zero theorem this work extends is due to **L. Rodina**
(arXiv:2406.04234); the framing of the locality programme owes much
to ongoing conversations with Rodina (professor) and with **Giuseppe**
(tutor). Any errors in the satellites here are the author's own.
