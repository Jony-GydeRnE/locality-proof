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
| **Foundational lemma 1.** Self-contained proof of $B \equiv 0 \Rightarrow B_i \equiv 0$ on each cyclic 1-zero zone (homogeneity decomposition). | `paper/Details missing in Hidden zero -> unitarity Rodina proof/Proof of Rodina claim B=0->B_i = 0/subset_vanishing.pdf` and `notes/13. alternate-proof-more-general.pdf` |
| **Foundational lemma 2.** Self-contained proof that the d-subset uniqueness argument works for all $n$. | `paper/Details missing in Hidden zero -> unitarity Rodina proof/details for Rodina D subset argument/dsubset_uniqueness.pdf` and `notes/18. d subset rigorous proof.pdf` |
| **Locality at $n = 4, 5, 6$ by hand.** Cyclic 1-zeros force every non-local coefficient to vanish via the Step-1 kill mechanism — 100% of non-locals killed at $n \le 6$. | `paper/step1 Kill technique and statistics/locality_unitarity_v5.pdf` |
| **Step-1 kill statistics across $n$.** 100% at $n \le 6$, 99.7–99.9% at $n = 7, 8, 9$. | `paper/step1 Kill technique and statistics/step-one-kill-statistics/step1_statistics.pdf` |
| **Locality at $n = 7, 8$.** Every Step-1 survivor dies via a depth-1 Laurent cascade — verified analytically (representative cases) and computationally (complete enumeration). | `paper/step2 Laurent series for hard kills/cascade_n7.pdf` (analytic write-up); `computations/step3_laurent/cascade_n7/` ($n=7$ enumeration); `computations/step4_laurent_block_analysis/n8/` ($n=8$ enumeration) |
| **Full nullspace dim $= 1$ at $n = 5, 6, 7$.** End-to-end SVD verification of the complete constraint system. | `computations/full_nullspace_verification/` |
| **Locality at $n = 9$ (new).** Step-1 kills 904,752 of 905,763 non-locals directly. The 1,011 survivors form 113 cyclic orbits — 23 close by independent depth-1 cascades, and the remaining **90 form a coupled cluster whose depth-1 fingerprint matrix has rank 90 (full)**, forcing every cluster coefficient to vanish jointly. First $n$ at which a block-rule mechanism is required. | `paper/step3 block rule/step3_block_rule_n9_REVISED.pdf` (write-up) and `paper/step3 block rule/step3_block_rule_README_REVISED.md`; per-orbit consolidated audit (all 113 reps) at `computations/step4_laurent_block_analysis/n9/outputs/n9_locality_status.md` |
| **Structural lemma.** Step-1 *never* removes a triangulation at any zone. Two-page exhaustive proof. | `paper/step1 Kill technique and statistics/step-one-doesnt-kill-triangulations/triangulation_layer0.pdf` |
| **All-$n$ working conjecture.** Locality follows from a rank condition on the depth-1 fingerprint matrix on each connected component of the cousin graph (verified through $n = 9$). | §6 below; refined statement in `paper/step3 block rule/step3_block_rule_README_REVISED.md` |
| **Reproducibility.** Every script that produced a number in this repo is committed alongside its raw output and a folder README. | `computations/` (top-level guide); each subfolder has its own README |
| **Working notes (~70 PDFs, condensing into a ~30-page undergrad-accessible companion).** | `notes/` — load-bearing: notes 13–17; the d-subset proof is note 18 |

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
The lemma in `paper/step1 Kill technique and statistics/step-one-doesnt-kill-triangulations/` proves the
strictly stronger statement that *no triangulation $T$ is in
$\mathcal K_{r,r+2}$ at any zone $r$* — ruling out both the originally
conjectured $T \subseteq F_r$ case and the case where $T$ contains a
substitute without its companion. The proof is two pages, three
exhaustive cases, using only the polygon edge $(r, r{+}1)$.

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
> alone.** Unitarity (all triangulation coefficients equal one common
> scalar) is handled separately by the d-subset paper in
> `paper/Details missing in Hidden zero -> unitarity Rodina proof/details for Rodina D subset argument/`.

At $n = 9$ Step 1 alone already kills **99.888%** of all non-triangulation
multisets. The hope is to prove the survivor fraction drops to $0$ fast as
$n$ grows, and that the residual "hard kills" can be handled uniformly by
the Laurent cascade.

---

## 4. How to read this repo (the guided tour)

Follow this order. Each section ends with a pointer to the relevant file
or subfolder. General notes are found in the Notes. This content takes one from trace phi 3 lagrangian and Feynman diagrams as triangulations to the the detailed proof that unitarity emerges at all n via hidden zeroes, following Rodina's paper https://arxiv.org/abs/2406.04234.  

### 4.1 Start with the geometric story
For readers who have never seen scattering amplitudes:

→ `paper/elementary-geometric-background-and-understanding/geometric_story.pdf`

The gentlest possible introduction to the polygon ↔ Feynman-diagram
dictionary, what a "hidden zero" is geometrically, and a closed-timelike-curve
heuristic for why locality has to hold. No prior knowledge assumed.

### 4.2 Read the by-hand worked example (the current "core" placeholder)
Once the setup makes sense:

→ `paper/step1 Kill technique and statistics/locality_unitarity_v5.pdf`

The full kill mechanism worked through at $n = 4, 5, 6$ explicitly, so you
can see Step 1 and Step 2 in action. This is the placeholder while the
core paper is still being written.

### 4.3 Read the foundational lemmas
The kill mechanism rests on two foundational lemmas:

→ `paper/Details missing in Hidden zero -> unitarity Rodina proof/Proof of Rodina claim B=0->B_i = 0/subset_vanishing.pdf`
   (homogeneity decomposition: $B|_{\mathcal Z_r}=0 \Rightarrow B_i|_{\mathcal Z_r}=0$)

→ `paper/Details missing in Hidden zero -> unitarity Rodina proof/details for Rodina D subset argument/dsubset_uniqueness.pdf`
   (uniqueness of the master substitution from the d-subset structure)

Together these justify applying the kill move term-by-term.

### 4.4 Read the Laurent cascade for hard kills
For why Step 1 alone is not enough at $n \ge 7$:

→ `paper/step2 Laurent series for hard kills/cascade_n7.pdf`

Proves all 7 of the $n=7$ "fish" die via a depth-1 Laurent cascade.
At $n = 9$ singleton cascade fails for some orbits and the **block-rule**
generalisation kicks in:

→ `paper/step3 block rule/step3_block_rule_n9_REVISED.pdf`

The 90-orbit cluster around orbit 22 has a depth-1 fingerprint matrix
of full rank — collectively kills every cluster orbit including the
4 singletons-cascade failures.

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
│   ├── Details missing in Hidden zero -> unitarity Rodina proof/
│   │   ├── Proof of Rodina claim B=0->B_i = 0/  §3 — homogeneity decomposition lemma
│   │   └── details for Rodina D subset argument/ §3 — d-subset uniqueness lemma
│   ├── step1 Kill technique and statistics/   §4 + §5 — Step-1 layer-0 kill
│   │   ├── locality_unitarity_v5               worked examples (n=4,5,6)
│   │   ├── step-one-doesnt-kill-triangulations local-survival lemma (Item 2 closed)
│   │   └── step-one-kill-statistics            global kill rates + asymptotic
│   ├── step2 Laurent series for hard kills/   §6 — singleton Laurent cascade (n=7)
│   └── step3 block rule/                      §7 — block-rule cluster kill (n=9)
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
| 1 | Foundational lemma $B = 0 \Rightarrow B_i = 0$ | draft in `paper/Details missing.../Proof of Rodina claim B=0->B_i = 0/` |
| 2 | No triangulation is killed by Step 1 at any zone | **closed** — `paper/step1 Kill technique and statistics/step-one-doesnt-kill-triangulations/` |
| 3 | Use Item 2 to identify exactly which indices Step 1 kills | now mechanical (Item 2 closed) |
| 4 | Survivor fraction $\to 0$ fast as $n$ grows | data through $n=9$; asymptotic conjecture open |
| 5 | Numerical verification at $n = 7, 8, 9$ | done ($n=7,8,9$); see `computations/step3_laurent/` and `computations/step4_laurent_block_analysis/` |
| 6 | Dual experiment ($X_{13}$ never special) | done — `computations/step1_layer0_kill/dual_X13_never_special/` |
| 7 | Conditional theorem: # survivors $= C_{n-2}$ | reduces to Item 10 |
| 8 | Worked examples at $n = 5, 6$ | done — `paper/step1 Kill technique and statistics/locality_unitarity_v5.pdf` |
| 9 | Undergraduate guide derived from the notes | in progress; `geometric_story.pdf` is the seed |
| 10 | Per-orbit depth-1 Laurent cascade closes all Step-1 survivors | **refined**; singleton verified at $n=7,8$; at $n=9$ generalises to BLOCK-RULE (cluster rank = full on the 90-orbit cluster). **Locality at $n=9$ is fully proven** (Step-1 + 23 single-orbit cascades + block-rule cluster kill ⟹ every non-triangulation coefficient is zero). |
| 11 | Triangulation $c$-values collapse to a single common scalar at $n \ge 9$ (UNITARITY) | handled separately by the **d-subset argument** in `paper/Details missing.../details for Rodina D subset argument/` and `notes/18.` — independent of the cascade machinery. (An earlier flip-graph framing of this question turned out to be a misframing — the 4 "untouched" Step-2 components at $n=9$ are triangulation components, which are LOCAL terms that should survive Step-1 by design, not be killed.) |

The updated dependency chain is:

**Items 1 + 2 (closed) ⟹ Item 3 (mechanical) ⟹ Item 10 (locality, via block-rule cluster kill) + Item 11 (unitarity, via d-subset) ⟹ locality + unitarity theorem.**

---

## 7. Contributing

### The non-negotiable rule

**If you produce experimental results that inform the research
direction of this proof, the code that produced them AND the
result files MUST land in this repo, in a folder under
[`computations/`](computations/).** No load-bearing claim about
locality, kill rates, cascade depth, cluster rank, etc. is
allowed in the paper or in any working note unless a reader can
go to a folder under `computations/` and find:

- the script that produced it (well commented per the rules below),
- the raw output (`.txt` / `.json` / `.md`) the script wrote,
- a README explaining what was tested, how long it took, what command
  ran, what the result means mathematically, and what remains unknown.

This is so that a mathematician who can't read code, but knows the
proof, can still verify *what was actually computed* by reading the
README and the result files. And so that another contributor (or a
future you) can re-run any experiment without having to reconstruct
the setup.

### Code style for the scripts themselves

If you write code for this repo (verifier scripts, enumerators,
symbolic checks, figure generators), follow these conventions so
that mathematicians and physicists can read the code without prior
Python fluency, and so that programmers can read it without prior
knowledge of the proof:

1. **Every function has a docstring with two parts:**
   - *LOGIC* — what algorithm is implemented (in CS / Python terms).
   - *PHYSICS / MATHEMATICS* — what the algorithm computes in the
     language of the proof (chords, zones, kill mechanism, etc.).
   See `computations/step4_laurent_block_analysis/n8/scripts/cascade_kill_n8.py`
   for the canonical example.
2. **Every script has a top-of-file block comment** stating its purpose,
   the conjecture it tests, the inputs/outputs, and any conventions that
   differ from the rest of the repo. If conventions are shared with
   another file (e.g. chord normalization), say so explicitly.
3. **Add a `README.md` to every new folder.** State what the folder is
   for, what files it contains, how to run them, and what results to
   expect. Headline numbers (kill rates, survivor counts, etc.) belong
   in the README, not buried in the script's stdout.
4. **No silent constants.** Every magic number (Laurent order to expand
   to, max number of survivors to check, etc.) is named and commented.
5. **If the script prints results, also write them to a text file** so
   they are reviewable without re-running.
6. **Push the code and the results together.** Don't push a paper claim
   like "all 100 survivors die" without simultaneously pushing the
   verifier and its trace; don't push a verifier without its
   results; don't put either of them outside `computations/`.

The goal is that any working physicist or mathematician who knows the
proof but has never seen this codebase should be able to open any script
and follow the logic in one read, and any contributor adding a new
experiment should know exactly what conventions to match.

---

## Acknowledgements

The hidden-zero theorem this work extends is due to **L. Rodina**
(arXiv:2406.04234); the framing of the locality programme owes much
to ongoing conversations with Rodina (professor) and with **Giuseppe**
(tutor). Any errors in the satellites here are the author's own.
