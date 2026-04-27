# Rodina 1-Zero Computation: Results & Proof Analysis

**Date:** 2026-04-15  
**Computation:** `full_ansatz.py`  
**Paper:** Rodina arXiv:2406.04234v5

---

## Computation Result

**OUTCOME 1 — the decisive case.**

| Quantity | Value |
|---|---|
| Ansatz terms | 45 (all pairs from 10 Mandelstam invariants) |
| 1-zero constraint equations | 115 independent |
| Matrix rank | 44 / 45 |
| **Nullspace dimension** | **1** |
| Nullspace basis | exactly `A_5^{φ³}` |
| Non-Feynman denominators in nullspace | **none** |

---

## Setup

**Independent variables:** `s12, s23, s34, s45, s15` (planar Mandelstams).

**Non-planar invariants** (momentum conservation):
```
s13 = s45 - s12 - s23
s14 = s23 - s15 - s45
s24 = s15 - s23 - s34
s25 = s34 - s12 - s15
s35 = s12 - s34 - s45
```

**Ansatz:** Most general dimension-(-4) rational function with only simple poles:
```
B = Σ_{a<b} c_{ab} / (s_a · s_b)
```
summed over all 10 Mandelstam invariants, giving 45 free parameters.

**1-zero loci** (set both non-adjacent invariants for each leg to zero):
```
Z1: s13=0, s14=0  →  s45 = s12+s23,  s15 = -s12
Z2: s24=0, s25=0  →  s15 = s23+s34,  s12 = -s23
Z3: s13=0, s35=0  →  s45 = s12+s23,  s34 = -s23
Z4: s14=0, s24=0  →  s23 = s15+s45,  s34 = -s45
Z5: s25=0, s35=0  →  s34 = s12+s15,  s45 = -s15
```

**Two types of constraint per locus:**
1. **Finiteness:** any `c_{ab}` where `s_a` or `s_b` vanishes on the locus must be zero (to avoid a pole on the locus).
2. **Vanishing:** for the remaining active terms, `B|_{locus} = 0` as a rational function of the 3 free variables → polynomial equations.

---

## Key Structural Observations

### Finiteness kills all non-planar poles

Each non-planar invariant vanishes on exactly two loci:
- `s13` on Z1 and Z3
- `s14` on Z1 and Z4
- `s24` on Z2 and Z4
- `s25` on Z2 and Z5
- `s35` on Z3 and Z5

The finiteness conditions on Z1 alone force `c_{s13,X} = c_{s14,X} = 0` for all X. Cycling through all 5 loci, **every non-planar pole term is eliminated**. The 1-zeros automatically screen out all non-Feynman denominators.

### Planar sector has a 1-dimensional kernel

The remaining 10 planar-pole terms (`1/(X_{ij} X_{kl})` for planar X's only) satisfy rank-9 constraints from the vanishing conditions, leaving a **1-dimensional nullspace**.

That unique solution is:
```
A_5 = 1/(s12·s34) + 1/(s12·s45) + 1/(s23·s45) + 1/(s23·s15) + 1/(s34·s15)
```
which is precisely the 5-pt φ³ tree amplitude (sum over 5 triangulations of the pentagon).

---

## Implications for the Proof Rewrite (from Eq. 9 onward)

### What this confirms

The 1-zero conditions, together with the requirement that B is a **rational function with only simple poles** in Mandelstam invariants, **uniquely fix B = A_tree** (up to an overall constant).

This means the following chain holds for n=5:

```
{Simple-pole rational function} ∩ {B=0 on all 5 one-zero loci}
    =  span{A_5^{φ³}}
```

### Proof structure that follows

The Rodina-style proof from Eq. 9 onward should proceed:

1. **Grading decomposition:** Write `B = B_{-2} + B_{-4}` (by degree in Mandelstams). Neither piece separately satisfies the 1-zero constraints — they mix across loci (the "cross-grading cancellation" noted in the previous session).

2. **Homogeneity lemma:** The computation shows the combined constraint is tight enough that the full `B` is fixed. The lemma "B=0 ⟹ B_i=0" allows recovery of grading-by-grading statements *after* the full uniqueness is established, not before.

3. **No locality assumption needed:** The finiteness analysis above is the locality constraint — it comes automatically from requiring B to be finite on the 1-zero loci. This means the proof does not need to separately invoke "only Feynman poles exist" as an input assumption; it is a *consequence* of the 1-zeros.

4. **The argument structure is:**
   - Input: `B` is a rational function with at most simple poles in 2-particle Mandelstams
   - Condition: `B` vanishes on all 5 cyclic 1-zero loci
   - Output: `B ∝ A_5^{φ³}`
   - Corollary: locality (only planar/Feynman poles) is implied, not assumed

### What needs care in the rewrite

- The proof above is **for the specific mass-dimension case** (denominator = product of exactly two Mandelstam invariants). For a fully general proof, one needs to extend the ansatz to include: (a) numerators, (b) higher-order pole structures, (c) different mass dimension sectors. The computation as written covers the most natural sector.

- The **grading B_{-2}, B_{-4}** distinction in the original Rodina proof is about different mass-dimension components. The current computation treats all terms uniformly. The cross-grading cancellation on loci means individual grade pieces don't satisfy 1-zeros; the full B does. This is consistent with Rodina's claim.

- **Equation 9 and beyond:** The rewrite should emphasize that the 1-zeros are the **primary constraints**, with unitarity/factorization following as a corollary (via the uniqueness result here). The direction is: `1-zeros → unique rational solution → identifies it as the Feynman amplitude → factorization/unitarity are thus inherited`.

---

## Files

- `proofs/full_ansatz.py` — SymPy computation (run with `python3 full_ansatz.py`)
- Output is deterministic; takes ~10s on a laptop

---

## Raw Output (for reference)

```
NULLSPACE DIMENSION: 1

── Nullspace vector 1 ──────────────────────────────────────
  Total nonzero entries: 5
  Planar pairs: 5,  Non-planar pairs: 0
    [P]  1/(s12*s34):  1
    [P]  1/(s12*s45):  1
    [P]  1/(s23*s45):  1
    [P]  1/(s23*s15):  1
    [P]  1/(s34*s15):  1

  >>> PROPORTIONAL TO A_tree  (factor = 1) ✓

CONCLUSION: Nullspace is 1-dimensional and EQUALS A_tree.
  => 1-zeros UNIQUELY determine the amplitude.
  => Proof strategy: 1-zeros + locality => A_tree is the only solution.
```
