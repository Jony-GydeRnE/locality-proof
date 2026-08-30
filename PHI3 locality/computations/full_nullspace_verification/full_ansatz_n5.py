#!/usr/bin/env python3
"""
full_ansatz_n5.py — 5-pt phi^3 theory: 1-zero nullspace computation.

PURPOSE
-------
Determines whether the 5-pt generic mass-dimension -4 ansatz (all 45 pairs of
the 10 Mandelstam invariants) has a 1-dimensional nullspace under the five cyclic
1-zero constraints, and whether that unique solution equals the phi^3 tree amplitude.

This script computationally verifies Theorem 1 (Part I) of the paper:
  "A purely algebraic proof of Rodina's hidden-zero theorem"

KINEMATICS
----------
n = 5 massless particles with momenta p_1,...,p_5 satisfying
  p_i^2 = 0   (massless)
  sum_i p_i = 0   (momentum conservation)

We use 2-particle Mandelstam invariants:
  s_ij = (p_i + p_j)^2

Under masslessness + momentum conservation the 10 invariants s_ij (1<=i<j<=5)
are not all independent. We choose 5 planar variables as independent:
  s12, s23, s34, s45, s15

The remaining 5 are expressed via momentum conservation:
  s13 = s45 - s12 - s23
  s14 = s23 - s15 - s45
  s24 = s15 - s23 - s34
  s25 = s34 - s12 - s15
  s35 = s12 - s34 - s45

1-ZERO LOCI
-----------
For each external leg k, the cyclic 1-zero locus Z_k sets to zero the two
non-adjacent 2-particle invariants involving leg k:

  Z1: s13 = 0,  s14 = 0   (non-adjacent to legs 2,5 which are adjacent to 1)
  Z2: s24 = 0,  s25 = 0
  Z3: s13 = 0,  s35 = 0
  Z4: s14 = 0,  s24 = 0
  Z5: s25 = 0,  s35 = 0

Each locus is imposed as a substitution in the planar variables.

ANSATZ
------
The most general rational function with simple poles in 2-particle Mandelstams
and total mass dimension -4 is:
  B = sum_{a < b} c_{ab} / (s_a * s_b)
where the sum runs over all C(10, 2) = 45 pairs, giving 45 free parameters.

ALGORITHM
---------
For each locus Z_k:
  1. Finiteness: any term c_{ab}/(s_a * s_b) where s_a or s_b vanishes on Z_k
     would be singular; force c_{ab} = 0.
  2. Vanishing: the remaining finite terms must sum to zero as a polynomial
     identity in the 3 free variables on Z_k; extract all coefficient equations.

The resulting linear system in the c's is solved for its nullspace.

EXPECTED OUTPUT
---------------
- Matrix rank 44/45, nullspace dimension 1
- Unique nullspace vector = A_tree = 1/(s12*s34) + 1/(s12*s45)
                                    + 1/(s23*s45) + 1/(s23*s15) + 1/(s34*s15)
- Zero non-Feynman (non-planar) denominators in the nullspace
- Conclusion: 1-zeros uniquely determine the amplitude (up to normalization)

HOW TO RUN
----------
  python3 full_ansatz_n5.py

Requires: sympy (pip install sympy)
Approximate runtime: ~10-15 seconds on a laptop.

REFERENCES
----------
Rodina, arXiv:2406.04234v5
"""

from sympy import (symbols, expand, together, fraction, Poly,
                   Matrix, linear_eq_to_matrix, zeros, Mul, Integer,
                   simplify, Rational)
from itertools import combinations

print("=" * 65)
print("5-pt phi^3  1-zero nullspace computation")
print("=" * 65)

# ── Independent planar variables ──────────────────────────────────────────────
# These 5 are algebraically independent at generic kinematics.
s12, s23, s34, s45, s15 = symbols('s12 s23 s34 s45 s15')

# ── All 10 Mandelstam invariants ──────────────────────────────────────────────
# Non-planar invariants expressed via momentum conservation:
#   sum of all s_ij = 0  (from (sum p_i)^2 = 0 with p_i^2 = 0)
s13_expr = s45 - s12 - s23
s14_expr = s23 - s15 - s45
s24_expr = s15 - s23 - s34
s25_expr = s34 - s12 - s15
s35_expr = s12 - s34 - s45

labels = ['s12', 's23', 's34', 's45', 's15', 's13', 's14', 's24', 's25', 's35']
exprs  = [ s12,   s23,   s34,   s45,   s15,
           s13_expr, s14_expr, s24_expr, s25_expr, s35_expr ]
N = 10
planar_idx = list(range(5))   # indices 0..4 are the 5 planar invariants

print("\nMomentum conservation check:")
total = sum(exprs)
print(f"  sum(all s_ij) = {expand(total)}   (should be 0)")

# ── The known 5-pt phi^3 tree amplitude ───────────────────────────────────────
# A_5 = sum over the 5 non-crossing triangulations of the convex pentagon:
#   each term 1/(s_{ij} * s_{kl}) with {ij},{kl} non-crossing pairs
A_tree_symbolic = (1/(s12*s34) + 1/(s12*s45) + 1/(s23*s45)
                   + 1/(s23*s15) + 1/(s34*s15))

# ── Ansatz: all C(10, 2) = 45 terms with simple poles ─────────────────────────
pairs = list(combinations(range(N), 2))   # list of (i, j) index pairs with i < j
M = len(pairs)
print(f"\nAnsatz: {M} terms  (all C(10,2) pairs from {N} Mandelstam invariants)")
print("  B = sum_{a<b} c_{ab} / (s_a * s_b)")

c = symbols(f'c0:{M}')   # 45 free parameters c0...c44

# ── 1-zero loci ───────────────────────────────────────────────────────────────
# Each locus is a substitution that sets two non-planar invariants to zero.
# The substitution is derived by solving the two zero-conditions for two planar vars.
#
# Z1: s13=0, s14=0  =>  s45 = s12+s23,  s15 = -s12
# Z2: s24=0, s25=0  =>  s15 = s23+s34,  s12 = -s23
# Z3: s13=0, s35=0  =>  s45 = s12+s23,  s34 = -s23
# Z4: s14=0, s24=0  =>  s23 = s15+s45,  s34 = -s45
# Z5: s25=0, s35=0  =>  s34 = s12+s15,  s45 = -s15
loci = [
    ({s45: s12+s23, s15: -s12},    ['s13', 's14']),  # Z1
    ({s15: s23+s34, s12: -s23},    ['s24', 's25']),  # Z2
    ({s45: s12+s23, s34: -s23},    ['s13', 's35']),  # Z3
    ({s23: s15+s45, s34: -s45},    ['s14', 's24']),  # Z4
    ({s34: s12+s15, s45: -s15},    ['s25', 's35']),  # Z5
]

# ── Verify locus derivations ──────────────────────────────────────────────────
print("\nLocus verification (vanishing invariants should evaluate to 0):")
for k, (sub, zero_lbls) in enumerate(loci):
    for lbl in zero_lbls:
        val = exprs[labels.index(lbl)].subs(sub)
        status = "✓" if expand(val) == 0 else "✗ ERROR"
        print(f"  Z{k+1}: {lbl} = {expand(val)}   {status}")

# ── Verify A_tree vanishes on all loci ────────────────────────────────────────
print("\nVerifying known A_5 vanishes on each 1-zero locus:")
for k, (sub, zero_lbls) in enumerate(loci):
    val = simplify(A_tree_symbolic.subs(sub))
    status = "✓ zero" if val == 0 else "✗ NONZERO"
    print(f"  Z{k+1}:  A_5 = {val}   {status}")

# ── Build linear constraints ──────────────────────────────────────────────────
print("\nBuilding linear constraints from 1-zero conditions...")
all_eqs = []   # list of linear expressions in c that must equal 0

for k, (sub, zero_lbls) in enumerate(loci):
    zero_idx = [labels.index(z) for z in zero_lbls]

    # CONSTRAINT TYPE 1: FINITENESS
    # Any term c_{ab}/(s_a * s_b) where s_a or s_b vanishes on Z_k would be
    # singular on the locus (blow up to infinity). Since B must be finite and
    # actually zero on Z_k, the coefficient of such a term must be zero.
    fin_count = 0
    for t, (i, j) in enumerate(pairs):
        if i in zero_idx or j in zero_idx:
            all_eqs.append(c[t])
            fin_count += 1

    # CONSTRAINT TYPE 2: VANISHING ON LOCUS
    # The remaining (non-singular) terms must sum to zero as a rational
    # function identity in the 3 free variables. Equivalently, the numerator
    # of the combined fraction must be the zero polynomial.
    sub_exprs = [e.subs(sub) for e in exprs]
    active = [(t, sub_exprs[i] * sub_exprs[j])
              for t, (i, j) in enumerate(pairs)
              if i not in zero_idx and j not in zero_idx]

    # The 3 free variables on this locus (those not fixed by the substitution)
    free = [v for v in [s12, s23, s34, s45, s15] if v not in sub]

    if active:
        # Combine: B|_{locus} = sum c[t] / denominator(t)
        # Put over common denominator and extract numerator polynomial
        B_locus = sum(c[t] / d for t, d in active)
        num, _den = fraction(together(B_locus))
        num_exp = expand(num)

        # Extract polynomial coefficients (each must vanish independently)
        try:
            poly = Poly(num_exp, *free, domain='EX')
            coeffs_list = poly.coeffs()
        except Exception:
            coeffs_list = []
            for mono in Poly(num_exp, *free).monoms():
                term = num_exp
                for var, exp in zip(free, mono):
                    term = term.coeff(var, exp)
                coeffs_list.append(expand(term))

        van_count = 0
        for coef in coeffs_list:
            coef_e = expand(coef)
            if coef_e != 0:
                all_eqs.append(coef_e)
                van_count += 1

        print(f"  Z{k+1}: {fin_count} finiteness eqs,  {van_count} vanishing eqs,  "
              f"{len(active)} active terms,  free vars: {[str(v) for v in free]}")
    else:
        print(f"  Z{k+1}: {fin_count} finiteness eqs,  0 active terms")

# ── Deduplicate equations ─────────────────────────────────────────────────────
seen = set()
unique_eqs = []
for eq in all_eqs:
    s = str(eq)
    if s not in seen and eq != 0:
        seen.add(s)
        unique_eqs.append(eq)

print(f"\nTotal equations (before dedup): {len(all_eqs)}")
print(f"Unique non-trivial equations:  {len(unique_eqs)}")

# ── Build and solve the linear system ─────────────────────────────────────────
print("\nBuilding coefficient matrix...")
A_mat, _ = linear_eq_to_matrix(unique_eqs, list(c))
rk = A_mat.rank()
print(f"Matrix shape: {A_mat.shape},  rank = {rk}")
print(f"Nullspace dimension (M - rank): {M - rk}")

# ── Compute nullspace ─────────────────────────────────────────────────────────
print("\nComputing nullspace (may take a moment)...")
ns = A_mat.nullspace()
print(f"\n{'='*65}")
print(f"NULLSPACE DIMENSION: {len(ns)}")
print(f"{'='*65}")

# ── Build the A_tree coefficient vector for comparison ────────────────────────
# A_5 = 1/(s12*s34) + 1/(s12*s45) + 1/(s23*s45) + 1/(s23*s15) + 1/(s34*s15)
A_tree_pairs_lbls = [('s12','s34'), ('s12','s45'), ('s23','s45'),
                     ('s23','s15'), ('s34','s15')]
A_tree_vec = [Integer(0)] * M
for la, lb in A_tree_pairs_lbls:
    ia, ib = sorted([labels.index(la), labels.index(lb)])
    t = pairs.index((ia, ib))
    A_tree_vec[t] = Integer(1)
A_tree_vec = Matrix(A_tree_vec)

print("\nKnown A_tree nonzero terms (these are the 5 non-crossing pentagon triangulations):")
for t in range(M):
    if A_tree_vec[t] != 0:
        i, j = pairs[t]
        print(f"  1/({labels[i]}*{labels[j]})")

# ── Analyse each nullspace vector ─────────────────────────────────────────────
print()
for vi, v in enumerate(ns):
    print(f"\n── Nullspace vector {vi+1} ──────────────────────────────────────")
    nonzero_t = [(t, v[t]) for t in range(M) if v[t] != 0]

    # Classify each nonzero term as planar (Feynman) or non-planar (non-Feynman)
    planar_terms = [(t, val) for t, val in nonzero_t
                    if pairs[t][0] in planar_idx and pairs[t][1] in planar_idx]
    nonplanar_terms = [(t, val) for t, val in nonzero_t
                       if not (pairs[t][0] in planar_idx and pairs[t][1] in planar_idx)]

    print(f"  Total nonzero entries: {len(nonzero_t)}")
    print(f"  Planar (Feynman) pairs: {len(planar_terms)},  Non-planar pairs: {len(nonplanar_terms)}")

    for t, val in nonzero_t:
        i, j = pairs[t]
        kind = "[P]" if (i in planar_idx and j in planar_idx) else "[NP]"
        print(f"    {kind}  1/({labels[i]}*{labels[j]}):  {val}")

    # Check proportionality to A_tree
    prop = None
    for t, val in planar_terms:
        if A_tree_vec[t] != 0:
            prop = simplify(val / A_tree_vec[t])
            break

    if prop is not None:
        residual = Matrix([simplify(v[t] - prop * A_tree_vec[t]) for t in range(M)])
        if all(r == 0 for r in residual):
            print(f"\n  >>> PROPORTIONAL TO A_tree  (factor = {prop}) ✓")
        else:
            nonzero_res = [(t, residual[t]) for t in range(M) if residual[t] != 0]
            print(f"\n  >>> NOT proportional to A_tree")
            print(f"      Residual nonzero at:")
            for t, r in nonzero_res:
                i, j = pairs[t]
                kind = "[P]" if (i in planar_idx and j in planar_idx) else "[NP]"
                print(f"        {kind}  1/({labels[i]}*{labels[j]}):  {r}")
    else:
        print(f"\n  >>> No overlap with A_tree (entirely non-Feynman?)")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print("SUMMARY")
print(f"{'='*65}")
print(f"Ansatz:          {M} terms (all C(10,2) pairs of Mandelstam invariants)")
print(f"1-zero loci:     5 cyclic legs, each zeroing 2 non-adjacent invariants")
print(f"Matrix rank:     {rk}  (out of {M} params)")
print(f"Nullspace dim:   {len(ns)}")
if len(ns) == 1:
    v = ns[0]
    prop = None
    for t in range(M):
        if A_tree_vec[t] != 0 and v[t] != 0:
            prop = simplify(v[t] / A_tree_vec[t])
            break
    residual = Matrix([simplify(v[t] - prop * A_tree_vec[t]) for t in range(M)])
    if all(r == 0 for r in residual):
        print(f"\nCONCLUSION: Nullspace is 1-dimensional and EQUALS A_tree.")
        print("  => 1-zeros UNIQUELY determine the amplitude (up to overall normalization).")
        print("  => This verifies Theorem 1 of the paper.")
    else:
        print(f"\nCONCLUSION: Nullspace is 1-dimensional but NOT A_tree.")
        print("  => Setup error — check loci or ansatz.")
elif len(ns) > 1:
    print(f"\nCONCLUSION: Nullspace has dim {len(ns)} — too large.")
    print("  Check if any extra vector involves non-Feynman denominators.")
else:
    print(f"\nCONCLUSION: Nullspace is EMPTY — something is wrong.")
    print("  A_tree should always be in the nullspace.")
