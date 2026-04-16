# Claude's Exploratory Proof Session — Hidden Zeros → Locality

## Context
This document records the step-by-step reasoning session between Jonathan Valenzuela and Claude that produced the algebraic proof of Rodina's hidden-zero theorem (Part I, n=5 base case). It serves as the "dirty work" behind the formal writeup in `paper/main.tex`.

## Key Results Established

### 1. Nullspace Computation
- 45-term generic ansatz (all mass-dim -2 rational functions in 10 Mandelstam invariants)
- 5 cyclic 1-zero loci imposed → 115 constraint equations
- Matrix rank: 44/45, nullspace dimension = 1
- Unique solution = φ³ tree amplitude A_tree
- Zero non-Feynman denominators survive

### 2. Two-Stage Locality Mechanism
- **Stage 1 (Finiteness):** On each locus Z_k, non-adjacent invariants vanish. Terms with those in the denominator blow up → coefficients forced to zero. Cycling all 5 loci kills all 35 non-planar terms.
- **Stage 2 (Vanishing):** Remaining 10 planar terms constrained by rank-9 system → 1 free direction = A_tree.

### 3. Independence of 1-Zeros
- Only 3 of 5 cyclic 1-zeros are independent (legs 1,2,3 suffice; legs 4,5 redundant)
- New observation not present in Rodina's original paper

### 4. The Non-Local Survivor
After legs 1+2: 2-parameter family c₁₀·A_tree + c₁₁·N where:
N = 1/(X₁₃·X₂₄) + 1/(X₁₃·X₂₅) + 1/(X₂₄·X₂₅) + 1/X₁₃²
Leg 3 forces c₁₁ = 0, killing the non-local piece. This is the precise mechanism by which the third 1-zero enforces locality.

### 5. Correct Logical Chain for Proof
Simple-pole rational B + B=0 on all 5 cyclic 1-zeros
  ⟹ locality (all non-Feynman poles vanish) [from finiteness]
  ⟹ B ∝ A_tree [from vanishing, uniqueness]
  ⟹ unitarity/factorization inherited

### 6. Cross-Grading Behavior
B₋₂ and B₋₁ do NOT satisfy 1-zeros individually — they cancel cross-grading on the locus. The homogeneity lemma (B=0 ⇒ Bᵢ=0) recovers grading-by-grading statements.

### 7. Leg-by-Leg Breakdown
- **Leg 1:** 9 equations, kills c₇,c₉,c₁₂,c₁₃,c₁₅ (crossing chords + squared propagators). 6 free params remain.
- **Leg 2:** 4 new equations, kills c₆,c₁₄, links c₅=c₁₀, c₈=c₁₀+c₁₁. 2 free params remain (c₁₀, c₁₁).
- **Leg 3:** 1 equation, forces c₁₁=0. 1 free param (c₁₀) = overall normalization.
- **Legs 4,5:** Redundant, contribute no new constraints.

## Implications for Part II (Induction)
The homogeneity lemma + dimensional reduction via 1-zeros should allow extending n=5 to general n. On Z_k, weight-i terms reduce to (n-1)-point ansatz. Induction hypothesis → restricted ansatz ∝ A_{n-1}^tree. BCFW recursion reconstructs A_n^tree.
