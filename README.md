# Locality and Unitarity of Tr(φ³) Tree Amplitudes from Hidden Zeros

> A working repository, by **Jonathan Valenzuela**, building toward an algebraic
> proof that the cyclic hidden 1-zero conditions uniquely determine the
> Tr(φ³) tree amplitudes — locality and unitarity emerge as *consequences*,
> not assumptions.

This README is the entry point. Read it top to bottom: it states the problem
precisely, gives the proof strategy in plain words, summarizes the current
status, and walks through every folder in the repo and how it supports the
conjecture.

---

## 1. The problem we want to prove

Label the legs of an n-point amplitude cyclically by 1,…,n. Write the
planar Mandelstam invariants $X_{ij} = (p_i+\dots+p_{j-1})^2$ for the
$d = n(n-3)/2$ chords of the n-gon. Let

$$B \;=\; \sum_{w}\frac{a_w}{\prod_{c\in w} X_c}$$

be the most general rational function of mass dimension $-2(n-3)$ with
at most simple poles on planar invariants ($w$ runs over size-$(n-3)$
multisets of chord indices). Define the cyclic hidden 1-zero locus

$$\mathcal{Z}_r \;:=\; \bigl\{c_{r,r+2}=c_{r,r+3}=\dots=c_{r,r+n-2}=0\bigr\},$$

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

The **central open question** of the upcoming core paper is purely
combinatorial:

> **Conjecture (Item 2).** For each zone $\mathcal Z_r$, every
> $(n-3)$-element subset of $F_r$ contains a crossing pair of chords.
> Equivalently: no triangulation of the $n$-gon is contained in any $F_r$.

This conjecture protects the triangulation coefficients from being
accidentally killed by Step 1, and is the missing combinatorial piece
of the all-$n$ proof. If proven — and the numerical kill rate at $n=9$
(99.889%) suggests it is true and the survivor fraction goes to $0$ fast —
the rest of the proof falls into place via a counting argument
(see §6 below, Items 7 and 4).

---

## 3. Status by n

| $n$ | Step 1 alone | Step 1 + 2 | Step 1 + 2 + 3 | Full nullspace |
|---|---|---|---|---|
| 4 | trivial | ✓ | ✓ | dim 1 |
| 5 | all 10 non-locals killed | ✓ | ✓ | **dim 1 (verified)** |
| 6 | 129 / 151 killed | ✓ remaining 22 die by cyclic Step 2 | ✓ | **dim 1 (verified)** |
| 7 | 7 fish escape | 7 fish still escape | ✓ all die by Laurent cascade | **dim 1 (verified)** |
| 8 | 100 escape | … | conjectured, not yet checked | — |
| 9 | 1 011 escape (out of 906 192) | … | conjectured, not yet checked | — |

At $n = 9$ Step 1 alone already kills **99.889%** of all non-triangulation
multisets. The hope is to prove the survivor fraction drops to $0$ fast as
$n$ grows, and that the residual "hard kills" can be handled uniformly by
the Laurent cascade.

---

## 4. How to read this repo (the guided tour)

Follow this order. Each section ends with a pointer to the relevant file
or subfolder.

### 4.1 Start with the geometric story
For readers who have never seen scattering amplitudes:

→ `paper/elementary-geometric-background-and-understanding/geometric_story.pdf`

The gentlest possible introduction to the polygon ↔ Feynman-diagram
dictionary, what a "hidden zero" is geometrically, and a closed-timelike-curve
heuristic for why locality has to hold. No prior knowledge assumed.

### 4.2 Read the by-hand worked example (the current "core" placeholder)
Once the setup makes sense:

→ `paper/ex- proving locality&unitarity by hand/locality_unitarity_v5.tex`

The full kill mechanism worked through at $n = 4, 5, 6$ explicitly, so you
can see Step 1 and Step 2 in action. This is the placeholder while the
core paper is still being written.

### 4.3 Read the foundational lemma
The kill mechanism rests on a homogeneity decomposition:

→ `paper/B=0->B_i = 0/subset_vanishing.pdf`

This generalizes Rodina's eqs ~15–17: if $B$ vanishes on $\mathcal Z_r$,
then so does each homogeneous weight-component $B_i$ separately. Needed
to apply the kill move term-by-term.

### 4.4 Read the Laurent cascade for hard kills
For why Step 1 alone is not enough at $n \ge 7$:

→ `paper/Laurent series for hard kills/cascade_n7.pdf`

Proves all 7 of the $n=7$ "fish" die via a depth-1 Laurent cascade. The
companion script reproduces every step symbolically.

### 4.5 Run the experiments
For numerical verification, survivor enumeration, and figures:

→ `computations/` — see `computations/README.md`

The folders are numbered to match the proof's logical flow:
`step0_sanity` → `step1_layer0_kill` → `step2_equate` → `step3_laurent`,
plus end-to-end checks (`full_nullspace_verification`) and figures
(`survivor_gallery`).

### 4.6 Browse the chronological notes
The complete working notebook (~50 PDFs):

→ `notes/`

Notes 1–17 are the topical PDFs (#1 = Feynman diagrams & triangulations,
#17 = hidden zeros ↔ non-local interactions). Notes 13–17 are the most
load-bearing for the current direction; the rest are dated session notes
from May 2025 onwards.

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
├── paper/                                   ← in-progress write-ups (one per §)
│   ├── README.md                              guided tour of paper/
│   ├── elementary-geometric-background…/      §1 — intuition for newcomers
│   ├── B=0->B_i = 0/                          §3 — foundational lemma
│   ├── ex- proving locality&unitarity by hand/  §4 — worked examples (n=4,5,6)
│   │                                            (also the standing core-paper placeholder)
│   └── Laurent series for hard kills/         §6 — Step-3 cascade at n=7
│
├── computations/                            ← experiments backing each proof step
│   ├── README.md                              guided tour of computations/
│   ├── step0_sanity/                          §2 — A_n^tree vanishes on each Z_r
│   ├── step1_layer0_kill/                     §4 — Step-1 kill enumeration + dual
│   ├── step2_equate/                          §4 — Step-2 equivalence classes
│   ├── step3_laurent/                         §6 — Laurent cascade at n=7
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
| 1 | Foundational lemma $B = 0 \Rightarrow B_i = 0$ | draft in `paper/B=0->B_i = 0/` |
| 2 | Every $(n-3)$-subset of $F_r$ contains a crossing pair | conjectured, **central open** |
| 3 | Use Item 2 to identify exactly which indices Step 1 kills | depends on Item 2 |
| 4 | Survivor fraction $\to 0$ fast as $n$ grows | quantitative asymptotic |
| 5 | Numerical verification at $n = 7, 8, 9$ | done — see `computations/` |
| 6 | Dual experiment ($X_{13}$ never special) | done — see `computations/step1_layer0_kill/dual_X13_never_special/` |
| 7 | Conditional theorem: # survivors $= C_{n-2}$ | proof strategy clear once Item 2 holds |
| 8 | Worked examples at $n = 5, 6$ | done — see `paper/ex- proving locality&unitarity by hand/` |
| 9 | Undergraduate guide derived from the notes | in progress; `geometric_story.pdf` is the seed |

The dependency chain is: **Items 1 + 2 ⟹ Item 3 ⟹ Item 7 ⟹ locality + unitarity theorem.**

---

## Author

**Jonathan Valenzuela**, 2026. Contact: see commit metadata.
