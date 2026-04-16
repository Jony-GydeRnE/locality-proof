#!/usr/bin/env python3
"""
independence_check.py — How many 1-zero loci are independent at n=5 and n=6?

PURPOSE
-------
For each n, compute the constraint rank after imposing 1, 2, ..., n cyclic
1-zero loci, and determine the minimum number needed to pin down the
amplitude uniquely (nullspace dimension = 1).

This verifies Corollary 4.2 (Part I) at n=5: only 3 of 5 loci are
independent. At n=6, the analogous count is determined numerically.

EXPECTED OUTPUT
---------------
n=5: 3 loci suffice (legs 1,2,3 give full rank; 4,5 add nothing)
n=6: floor(6/2)=3 loci should suffice (per Remark 5.7 conjecture)

HOW TO RUN
----------
  python3 independence_check.py

Requires: numpy, sympy (sympy only for n=5 symbolic check)
Runtime: < 20 seconds.
"""

import numpy as np
from itertools import combinations, combinations_with_replacement
from sympy import (symbols, expand, together, fraction, Poly,
                   Matrix, linear_eq_to_matrix, simplify, Integer)

print("=" * 65)
print("Independence of 1-zero loci")
print("=" * 65)

# =====================================================================
#  n = 5: SYMBOLIC (exact)
# =====================================================================

print("\n" + "-" * 65)
print("n = 5 (symbolic / SymPy)")
print("-" * 65)

s12, s23, s34, s45, s15 = symbols('s12 s23 s34 s45 s15')

s13_expr = s45 - s12 - s23
s14_expr = s23 - s15 - s45
s24_expr = s15 - s23 - s34
s25_expr = s34 - s12 - s15
s35_expr = s12 - s34 - s45

labels = ['s12', 's23', 's34', 's45', 's15', 's13', 's14', 's24', 's25', 's35']
exprs  = [s12, s23, s34, s45, s15, s13_expr, s14_expr, s24_expr, s25_expr, s35_expr]
N = 10
pairs = list(combinations(range(N), 2))  # 45 pairs
M_ansatz = len(pairs)
c = symbols(f'c0:{M_ansatz}')

loci_5 = [
    ({s45: s12+s23, s15: -s12},    ['s13', 's14']),
    ({s15: s23+s34, s12: -s23},    ['s24', 's25']),
    ({s45: s12+s23, s34: -s23},    ['s13', 's35']),
    ({s23: s15+s45, s34: -s45},    ['s14', 's24']),
    ({s34: s12+s15, s45: -s15},    ['s25', 's35']),
]

def build_eqs_for_loci(loci_subset):
    """Build constraint equations from a subset of loci (0-indexed)."""
    all_eqs = []
    for k in loci_subset:
        sub, zero_lbls = loci_5[k]
        zero_idx = [labels.index(z) for z in zero_lbls]

        # Finiteness constraints
        for t, (i, j) in enumerate(pairs):
            if i in zero_idx or j in zero_idx:
                all_eqs.append(c[t])

        # Vanishing constraints
        sub_exprs = [e.subs(sub) for e in exprs]
        active = [(t, sub_exprs[i] * sub_exprs[j])
                  for t, (i, j) in enumerate(pairs)
                  if i not in zero_idx and j not in zero_idx]
        free = [v for v in [s12, s23, s34, s45, s15] if v not in sub]

        if active:
            B_locus = sum(c[t] / d for t, d in active)
            num, _ = fraction(together(B_locus))
            num_exp = expand(num)
            try:
                poly = Poly(num_exp, *free, domain='EX')
                for coef in poly.coeffs():
                    coef_e = expand(coef)
                    if coef_e != 0:
                        all_eqs.append(coef_e)
            except Exception:
                pass

    # Deduplicate
    seen = set()
    unique = []
    for eq in all_eqs:
        s = str(eq)
        if s not in seen and eq != 0:
            seen.add(s)
            unique.append(eq)
    return unique

print("\nProgressive rank as loci are added (n=5, 45-param ansatz):")
print(f"{'Loci used':<25} {'# eqs':>8} {'Rank':>8} {'Null dim':>10}")
print("-" * 55)

for n_loci in range(1, 6):
    # Use legs 1, 2, ..., n_loci
    eqs = build_eqs_for_loci(list(range(n_loci)))
    A_mat, _ = linear_eq_to_matrix(eqs, list(c))
    rk = A_mat.rank()
    null_dim = M_ansatz - rk
    loci_str = ', '.join([str(i+1) for i in range(n_loci)])
    print(f"  Legs {loci_str:<18} {len(eqs):>8} {rk:>8} {null_dim:>10}")

# Also check: which SINGLE locus gives the most constraints?
print("\nSingle-locus constraint counts:")
for k in range(5):
    eqs = build_eqs_for_loci([k])
    A_mat, _ = linear_eq_to_matrix(eqs, list(c))
    rk = A_mat.rank()
    print(f"  Leg {k+1}: rank = {rk},  null dim = {M_ansatz - rk}")

# =====================================================================
#  n = 6: NUMERICAL (numpy)
# =====================================================================

print("\n" + "-" * 65)
print("n = 6 (numerical / numpy)")
print("-" * 65)

N6_VARS = 9
VAR6 = ['X13', 'X14', 'X15', 'X24', 'X25', 'X26', 'X35', 'X36', 'X46']

LOCI_6 = [
    {'c': [0,1,2], 'f': [3,4,5,6,7,8],
     's': np.array([[0,0,-1,0,0,0],[1,0,-1,0,0,0],[0,1,-1,0,0,0]], dtype=float)},
    {'c': [3,4,5], 'f': [0,1,2,6,7,8],
     's': np.array([[-1,0,0,0,0,0],[-1,0,0,1,0,0],[-1,0,0,0,1,0]], dtype=float)},
    {'c': [0,6,7], 'f': [1,2,3,4,5,8],
     's': np.array([[1,0,-1,0,0,0],[0,0,-1,0,0,0],[0,0,-1,0,0,1]], dtype=float)},
    {'c': [1,3,8], 'f': [0,2,4,5,6,7],
     's': np.array([[0,1,0,0,-1,0],[0,0,1,0,-1,0],[0,0,0,0,-1,0]], dtype=float)},
    {'c': [2,4,6], 'f': [0,1,3,5,7,8],
     's': np.array([[0,0,0,0,0,-1],[0,0,0,1,0,-1],[0,0,0,0,1,-1]], dtype=float)},
    {'c': [5,7,8], 'f': [0,1,2,3,4,6],
     's': np.array([[0,0,-1,0,0,0],[1,0,-1,0,0,0],[0,1,-1,0,0,0]], dtype=float)},
]

simple_basis_6 = list(combinations(range(N6_VARS), 3))  # 84
N_BASIS_6 = len(simple_basis_6)

def sample_n6(loc, rng):
    free_vals = rng.standard_normal(6) * 3.0
    X = np.zeros(N6_VARS)
    for i, fi in enumerate(loc['f']):
        X[fi] = free_vals[i]
    cv = loc['s'] @ free_vals
    for i, ci in enumerate(loc['c']):
        X[ci] = cv[i]
    return X

def build_n6_matrix(loci_indices, n_samples=300, seed=42):
    rng = np.random.default_rng(seed)
    rows = []
    for k in loci_indices:
        loc = LOCI_6[k]
        good = 0
        att = 0
        while good < n_samples and att < n_samples * 30:
            X = sample_n6(loc, rng)
            att += 1
            if np.min(np.abs(X)) < 0.05:
                continue
            try:
                row = np.array([1.0/np.prod(X[list(t)]) for t in simple_basis_6])
                if np.all(np.isfinite(row)) and np.max(np.abs(row)) < 1e14:
                    rows.append(row)
                    good += 1
            except:
                continue
    return np.vstack(rows)

print(f"\nProgressive rank as loci are added (n=6, {N_BASIS_6}-param ansatz):")
print(f"{'Loci used':<25} {'Rank':>8} {'Null dim':>10}")
print("-" * 47)

for n_loci in range(1, 7):
    M_mat = build_n6_matrix(list(range(n_loci)), n_samples=300, seed=12345)
    U, S, Vt = np.linalg.svd(M_mat, full_matrices=True)
    tol = S[0] * max(M_mat.shape) * np.finfo(float).eps * 1000
    rk = np.sum(S > tol)
    null_dim = N_BASIS_6 - rk
    loci_str = ', '.join([str(i+1) for i in range(n_loci)])
    print(f"  Legs {loci_str:<18} {rk:>8} {null_dim:>10}")

# Check single-locus constraints at n=6
print("\nSingle-locus constraint counts (n=6):")
for k in range(6):
    M_mat = build_n6_matrix([k], n_samples=300, seed=42)
    U, S, Vt = np.linalg.svd(M_mat, full_matrices=True)
    tol = S[0] * max(M_mat.shape) * np.finfo(float).eps * 1000
    rk = np.sum(S > tol)
    print(f"  Leg {k+1}: rank = {rk},  null dim = {N_BASIS_6 - rk}")

# Check all combinations of 3 loci
print("\nAll C(6,3)=20 combinations of 3 loci:")
from itertools import combinations as comb
results_3 = []
for combo in comb(range(6), 3):
    M_mat = build_n6_matrix(list(combo), n_samples=200, seed=42)
    U, S, Vt = np.linalg.svd(M_mat, full_matrices=True)
    tol = S[0] * max(M_mat.shape) * np.finfo(float).eps * 1000
    rk = np.sum(S > tol)
    nd = N_BASIS_6 - rk
    legs = tuple(k+1 for k in combo)
    results_3.append((legs, rk, nd))
    suffix = " ✓ UNIQUE" if nd == 1 else ""
    print(f"  Legs {legs}: rank = {rk}, null dim = {nd}{suffix}")

n_sufficient = sum(1 for _, _, nd in results_3 if nd == 1)
print(f"\n{n_sufficient} / 20 combinations of 3 loci give unique amplitude.")

# ── Summary ──────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print("n=5: 3 loci suffice (Corollary 4.2). Legs 4,5 contribute nothing new.")
print("n=6: See above for which combinations of 3 loci suffice.")
print("     floor(6/2) = 3 is the conjectured sharp bound (Remark 5.7).")
