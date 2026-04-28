# Hidden Zeros, Locality, and the Rodina Bootstrap

A reading guide and toolkit for understanding Rodina's hidden-zero uniqueness program (arXiv:2406.04234) and its connection to the locality emergence proof in Tr(φ³). Designed for someone new to the field who wants to come up to speed in a few weeks rather than a few months.

---

## What is this repo for?

Rodina's paper proves a striking claim: the hidden zeros discovered by Arkani-Hamed et al. uniquely determine the Tr(φ³) tree amplitude up to an overall normalization. From there, locality and tree-level unitarity follow as corollaries. The proof is short — about four pages of text — but very dense. Several arguments are stated only by example, and the connection between the algebraic mechanism and the underlying geometry is not always made explicit.

This repo collects:

1. **A guided path** through Rodina's paper, with the implicit steps filled in.
2. **Geometric intuition** for why the cyclic hidden zeros enforce locality, including the duality between Feynman diagrams and triangulations and a physical (causality-based) argument for why crossing chords cannot correspond to physical channels.
3. **By-hand proofs** of locality emergence at n = 4, 5, 6, with sympy/numerical verification at n = 6, 7, 8, 9.
4. **An expanded version** of Rodina's `B = 0 ⇒ B_i = 0` argument with minimal hypotheses and a one-line proof.
5. **A rigorous version** of Rodina's D-subset uniqueness argument (his eq. 22), generalized to arbitrary k and n.

The audience is a graduate student or advanced undergraduate with a working background in QFT and basic combinatorics, who wants to understand the bootstrap program in this corner of amplitudes well enough to contribute. The goal is to flatten the learning curve from months to one or two weeks.

---

## How to read this repo

### Path A — fastest route to understanding Rodina's main result

1. Start with **Rodina's paper** itself (`Hidden_zeros_are_equivalent...2406_04234v5.pdf`).
2. When you hit eqs. (15)–(18) and find the argument from "the only identities the X^(0) variables satisfy" terse, switch to **`subset_vanishing.pdf`**. This gives a clean, hypothesis-minimal proof of `B ≡ 0 ⇒ B_i ≡ 0` using only a single rescaling argument. After reading it, return to Rodina with a concrete picture of what eq. (15)–(18) is doing and why his Appendix A is needed only for the *enhanced*-scaling refinement, not for the basic equivalence.
3. When you hit eq. (22) and want to know how the "two cuts" argument generalizes beyond k = 2, switch to **`dsubset_uniqueness.pdf`**. This contains the explicit induction across the k+1 diagrams of a D-subset and the simultaneous-limit mechanism that isolates each adjacent pair of coefficients. **This PDF fills in details missing in Rodina's D-subset argument** — specifically, the argument that "the whole process can be carried out for each D-subset" is made precise, with an explicit proof that works for all k and n.

### Path B — for the locality emergence direction

1. Start with **`geometric_story.pdf`** for the polygon picture, the chord-birth interpretation of cyclic hidden zeros, the Feynman-diagram ↔ triangulation duality, and a four-part causality argument for why crossing chords must be killed.
2. Then read the by-hand proofs in **`locality_n7.pdf`** (despite the name, this contains the full n = 4, 5, 6 proofs) and **`localityunitarity.pdf`** (older draft — being rewritten; treat as a reference for the kinematic mesh setup).
3. Numerical verification at higher n: see `code/full_ansatz_n6.py`, `full_ansatz_n7.py`, and the result files `results_n6.txt`, `results_n7.txt`. For n = 8 and n = 9 see `nonlocal_n8_*`, `nonlocal_n9_*`.

### Path C — full background, slow read

The numbered notes `1__feynmandiagramsandtriangulations.pdf` through `18__d_subset_rigorous_proof.pdf` are a complete handwritten walkthrough, in pedagogical order, of every concept needed for Rodina's paper plus the locality program. They assume nothing beyond standard QFT.

---

## Files in this repo

### Polished PDFs (these are the entry points)

| File | What it is |
|---|---|
| `geometric_story.pdf` | Polygon, chord birth, Feynman ↔ triangulation duality, causality argument for locality. ~7 pages. |
| `subset_vanishing.pdf` | Clean, hypothesis-minimal proof of `B ≡ 0 ⇒ B_i ≡ 0`. Replaces Rodina's eqs. 15–18 with a single rescaling argument. ~2 pages. |
| `dsubset_uniqueness.pdf` | **Fills in details missing in Rodina's D-subset argument.** Generalizes Rodina's eq. 22 (worked only for k = 2) to a full induction over arbitrary k and n. ~3 pages. |

### Locality emergence proofs (by-hand)

| File | What it is |
|---|---|
| `locality_n7.pdf` | By-hand proofs of locality at n = 4, 5, 6. Despite the filename, n = 7 is treated only numerically. |
| `localityunitarity.pdf` | Older draft of the locality+unitarity paper. Contains the kinematic-mesh setup and Feynman-as-triangulation summary. **Being rewritten** — treat as reference, not as final form. |
| `5ptlocality.pdf`, `2ndorder5ptlocality.pdf`, `6ptlocality.pdf` | Step-by-step computations at fixed n. |
| `Locality_Emergence_Lemma` | Master formula and the local kill-test framing. |

### Numerical verification

| File | What it is |
|---|---|
| `full_ansatz_n6.py`, `full_ansatz_n7.py` | Numerical ansatz solvers at n = 6, 7. |
| `results_n6.txt`, `results_n7.txt` | Output of the above. |
| `nonlocal_n7_*`, `nonlocal_n8_*`, `nonlocal_n9_*` | Survivor lists from the kill-mechanism algorithm at higher n. |
| `step1_uncaught.py`, `step1_dual.py`, `step2_classes.py` | Helper scripts for the algorithm. |

### Numbered handwritten notes (complete pedagogical walkthrough)

| Notes | Topic |
|---|---|
| 1–4 | Feynman diagrams, Tr(φ³) Lagrangian, physical picture, triangulations |
| 5–8 | Walk-through of Rodina's preliminary sections, factorization, BCFW shifts, hidden-zero vs BCFW equivalence |
| 9 | Rodina's core theorem |
| 10–12 | Filling in details of Rodina's proof; novel revision of the BCFW-shift / hidden-zero correspondence; clean proof of `B = 0 ⇒ B_i = 0` |
| 13 | Alternate, more general proof of `B = 0 ⇒ B_i = 0` (the basis of `subset_vanishing.pdf`) |
| 14 | 6-point encroaching-locality proof (in progress) |
| 15 | Pattern structures in zeroes and chords |
| 16 | Cyclic shifts of zeroes |
| 17 | Hidden-zero / non-local interaction correspondence |
| 18 | **Rigorous D-subset proof** (the basis of `dsubset_uniqueness.pdf`) |

### Reference papers

| File | Paper |
|---|---|
| `Hidden_zeros_are_equivalent...2406_04234v5.pdf` | Rodina, the central paper |
| `Hidden_zeros_for_particlestring_amplitudes...2312_16282v3.pdf` | Arkani-Hamed et al., hidden zeros + δ-deformation |
| `geometricbackground1711_09102v2.pdf` | ABHY associahedron / kinematic-amplitude geometry |
| `LocalityEmerges45ptsintracephi3.pdf` | Earlier note on locality at low n |

---

## What's been filled in vs. Rodina

Three places in Rodina's paper are stated tersely or by example, and are made fully rigorous here:

1. **Eqs. 15–18** (`B ≡ 0 ⇒ B_m ≡ 0`, the "subset zero equivalence"). Rodina's argument runs through a structural representation `B = Σ c_ij {…}_ij` and an identification of `X^(∞)` and `X^(0)` relations (his Appendix A) before iterating from the top scaling. **`subset_vanishing.pdf`** replaces this with a one-paragraph proof using only a rescaling and the resulting Laurent polynomial identity. Appendix A is then needed only for the separate enhanced-scaling claim (`B_i ~ z^{i-1}` rather than `z^i`), not for the basic equivalence.

2. **Eq. 22** (D-subset uniqueness — all coefficients in a D-subset are equal). Rodina works the case k = 2 at n = 6 and asserts the general case follows analogously. **`dsubset_uniqueness.pdf`** gives the explicit induction across all k+1 diagrams of a D-subset, using two simultaneous limits per adjacent pair plus the master substitution `X_{2,j} = X_{1,j} - X_{1,3}`. The mechanism is uniform across k.

3. **Geometric intuition for locality**. Rodina's framework establishes locality as a corollary of uniqueness, but does not give a physical picture for *why* the hidden zeros must kill non-local terms. **`geometric_story.pdf`** supplies the polygon-with-chords picture, the master substitution as a chord-birth equation, the Feynman-diagram ↔ triangulation duality, and a four-part causality argument: a Feynman diagram with crossing chords would correspond to a closed time-oriented loop, hence to superluminal propagation, hence is forbidden.

---

## Status and open questions

- **Locality at general n:** the by-hand proofs at n = 4, 5, 6 plus numerics at n = 7, 8, 9 strongly suggest a general theorem. The combinatorial blocking-set lemma needed to close it is not yet proven. See `Locality_Emergence_Lemma` and notes 14, 17.
- **Rewrite of `localityunitarity.pdf`:** in progress. The new draft uses the master formula + local kill-test framing, with the n = 5, 6 cases as appendices and a stated blocking-set conjecture for general n.
- **Connection to the massive theory:** Rodina's recent work (arXiv:2601.16860, January 2026) shows that hidden zeros survive symmetry-controlled mass generation. The natural extension question is whether massive hidden zeros also uniquely fix amplitudes. This is the active thesis direction.

---

## Citation

If you use any of this material, please cite:

```
Rodina, L. "Hidden zeros are equivalent to enhanced ultraviolet scaling and lead to
unique amplitudes in Tr(φ³) theory." arXiv:2406.04234.
```

For the new proofs in this repo, a working citation form is:

```
[Author], "Geometric Picture and the D-Subset Uniqueness Theorem in Tr(φ³)."
Working notes, BIMSA, 2026.
```

---

*Last updated: April 2026. Comments and corrections welcome.*
