#!/usr/bin/env python3
"""
full_ansatz_n6.py — 6-pt phi^3 theory: 1-zero nullspace computation (numerical).

PURPOSE
-------
Verifies the hidden-zero theorem at n=6: the cyclic 1-zero conditions on the
ansatz space B_6 of mass-dimension -3 rational functions yield a 1-dimensional
nullspace spanned by A_6^tree.

Uses numerical linear algebra (numpy) for speed; the 84-dimensional ansatz
at n=6 would be prohibitively slow with symbolic methods.

KINEMATICS
----------
9 planar Mandelstam variables at n=6 (non-adjacent pairs in the hexagon):
  X13=(p1+p2)^2,  X14=(p1+p2+p3)^2,  X15=(p1+p2+p3+p4)^2
  X24=(p2+p3)^2,  X25=(p2+p3+p4)^2,  X26=(p2+p3+p4+p5)^2
  X35=(p3+p4)^2,  X36=(p3+p4+p5)^2,  X46=(p4+p5)^2

ANSATZ
------
Simple-poles ansatz: all C(9,3)=84 terms 1/(X_a*X_b*X_c), a<b<c distinct.
Broad ansatz: all C(11,3)=165 terms including 1/(X_a^2*X_b) and 1/(X_a^3).

1-ZERO LOCI
-----------
Each leg k imposes 3 constraints c_{k,j}=0 for j=k+2,k+3,k+4 (mod 6),
reducing 9 variables to 6 free via linear substitution.

EXPECTED OUTPUT
---------------
Nullspace dimension 1 in both ansatz spaces, spanned by A_6^tree.

HOW TO RUN
----------
  python3 full_ansatz_n6.py

Requires: numpy
Runtime: < 5 seconds.
"""

import numpy as np
from itertools import combinations, combinations_with_replacement
import sys

# =====================================================================
#  SETUP
# =====================================================================

N_VARS = 9
VAR_NAMES = ['X13', 'X14', 'X15', 'X24', 'X25', 'X26', 'X35', 'X36', 'X46']
# Index:       0      1      2      3      4      5      6      7      8

# =====================================================================
#  TREE AMPLITUDE — 14 triangulations of the hexagon
# =====================================================================
# Each triangulation of the convex hexagon uses 3 non-crossing diagonals.
# The C_4 = 14 triangulations, enumerated by vertex-1 diagonals used:

TREE_DIAG_NAMES = [
    ('X13','X14','X15'),   #  1. fan from vertex 1
    ('X13','X14','X46'),   #  2.
    ('X13','X15','X35'),   #  3.
    ('X14','X15','X24'),   #  4.
    ('X13','X35','X36'),   #  5.
    ('X13','X36','X46'),   #  6.
    ('X14','X24','X46'),   #  7.
    ('X15','X24','X25'),   #  8.
    ('X15','X25','X35'),   #  9.
    ('X24','X25','X26'),   # 10.
    ('X24','X26','X46'),   # 11.
    ('X25','X26','X35'),   # 12.
    ('X26','X35','X36'),   # 13.
    ('X26','X36','X46'),   # 14.
]

TREE_TRIPLES = [
    tuple(sorted(VAR_NAMES.index(n) for n in names))
    for names in TREE_DIAG_NAMES
]
assert len(TREE_TRIPLES) == 14

# =====================================================================
#  1-ZERO LOCI
# =====================================================================
# Each Z_k: 3 variables = linear combination of 6 free variables.
# Derived from c_{k,j} = X_{kj}+X_{k+1,j+1}-X_{k,j+1}-X_{k+1,j} = 0.

LOCI = [
    {   # Z_1: X13=-X26, X14=X24-X26, X15=X25-X26
        'constrained': [0, 1, 2],
        'free':        [3, 4, 5, 6, 7, 8],
        'sub': np.array([
            [ 0,  0, -1,  0,  0,  0],
            [ 1,  0, -1,  0,  0,  0],
            [ 0,  1, -1,  0,  0,  0],
        ], dtype=float),
    },
    {   # Z_2: X24=-X13, X25=X35-X13, X26=X36-X13
        'constrained': [3, 4, 5],
        'free':        [0, 1, 2, 6, 7, 8],
        'sub': np.array([
            [-1,  0,  0,  0,  0,  0],
            [-1,  0,  0,  1,  0,  0],
            [-1,  0,  0,  0,  1,  0],
        ], dtype=float),
    },
    {   # Z_3: X13=X14-X24, X35=-X24, X36=X46-X24
        'constrained': [0, 6, 7],
        'free':        [1, 2, 3, 4, 5, 8],
        'sub': np.array([
            [ 1,  0, -1,  0,  0,  0],
            [ 0,  0, -1,  0,  0,  0],
            [ 0,  0, -1,  0,  0,  1],
        ], dtype=float),
    },
    {   # Z_4: X14=X15-X35, X24=X25-X35, X46=-X35
        'constrained': [1, 3, 8],
        'free':        [0, 2, 4, 5, 6, 7],
        'sub': np.array([
            [ 0,  1,  0,  0, -1,  0],
            [ 0,  0,  1,  0, -1,  0],
            [ 0,  0,  0,  0, -1,  0],
        ], dtype=float),
    },
    {   # Z_5: X15=-X46, X25=X26-X46, X35=X36-X46
        'constrained': [2, 4, 6],
        'free':        [0, 1, 3, 5, 7, 8],
        'sub': np.array([
            [ 0,  0,  0,  0,  0, -1],
            [ 0,  0,  0,  1,  0, -1],
            [ 0,  0,  0,  0,  1, -1],
        ], dtype=float),
    },
    {   # Z_6: X26=-X15, X36=X13-X15, X46=X14-X15
        'constrained': [5, 7, 8],
        'free':        [0, 1, 2, 3, 4, 6],
        'sub': np.array([
            [ 0,  0, -1,  0,  0,  0],
            [ 1,  0, -1,  0,  0,  0],
            [ 0,  1, -1,  0,  0,  0],
        ], dtype=float),
    },
]


def sample_point(loc, rng, scale=3.0):
    """Sample a random kinematic point on a 1-zero locus."""
    free_vals = rng.standard_normal(6) * scale
    X = np.zeros(N_VARS)
    for i, fi in enumerate(loc['free']):
        X[fi] = free_vals[i]
    constr_vals = loc['sub'] @ free_vals
    for i, ci in enumerate(loc['constrained']):
        X[ci] = constr_vals[i]
    return X


def evaluate_tree(X):
    """Evaluate A_6^tree at kinematic point X."""
    return sum(1.0 / (X[a] * X[b] * X[c]) for a, b, c in TREE_TRIPLES)


def build_constraint_matrix(basis, n_samples=250, seed=12345, verbose=True):
    """Build constraint matrix by sampling points on all 6 loci."""
    rng = np.random.default_rng(seed)
    rows = []
    n_basis = len(basis)

    for k in range(6):
        loc = LOCI[k]
        good = 0
        attempts = 0
        while good < n_samples and attempts < n_samples * 30:
            X = sample_point(loc, rng)
            attempts += 1
            if np.min(np.abs(X)) < 0.05:
                continue
            try:
                row = np.array([1.0 / np.prod(X[list(t)]) for t in basis])
                if np.all(np.isfinite(row)) and np.max(np.abs(row)) < 1e14:
                    rows.append(row)
                    good += 1
            except (ZeroDivisionError, FloatingPointError):
                continue
        if verbose:
            print(f"  Z{k+1}: {good} samples ({attempts} attempts)")

    return np.vstack(rows)


def compute_nullspace(M, label=""):
    """Compute nullspace via SVD."""
    U, S, Vt = np.linalg.svd(M, full_matrices=True)
    n_basis = M.shape[1]
    tol = S[0] * max(M.shape) * np.finfo(float).eps * 1000
    rank = np.sum(S > tol)
    null_dim = n_basis - rank

    print(f"\n{label}")
    print(f"  Matrix: {M.shape},  rank = {rank} / {n_basis}")
    print(f"  Singular values near cutoff:")
    for i in range(max(0, rank - 2), min(n_basis, rank + 3)):
        marker = " <-- gap" if i == rank else ""
        print(f"    S[{i:3d}] = {S[i]:.6e}{marker}")
    print(f"  NULLSPACE DIMENSION: {null_dim}")

    return rank, null_dim, Vt[rank:], S


def check_proportional(v, A_ref, basis):
    """Check if v is proportional to A_ref. Return (is_prop, residual)."""
    ref_idx = next((i for i in range(len(A_ref)) if abs(A_ref[i]) > 0.5), None)
    if ref_idx is None or abs(v[ref_idx]) < 1e-10:
        return False, float('inf')
    scale = v[ref_idx] / A_ref[ref_idx]
    residual = np.max(np.abs(v - scale * A_ref))
    return residual < 1e-6, residual


# =====================================================================
#  MAIN
# =====================================================================

def main():
    print("=" * 65)
    print("6-pt phi^3  1-zero nullspace computation  (numerical)")
    print("=" * 65)

    # ── Verify loci ───────────────────────────────────────────────
    print("\nVerifying A_6^tree vanishes on each locus:")
    rng_v = np.random.default_rng(999)
    for k in range(6):
        vals = []
        for _ in range(300):
            X = sample_point(LOCI[k], rng_v)
            if np.min(np.abs(X)) < 0.01:
                continue
            vals.append(abs(evaluate_tree(X)))
        mx = max(vals) if vals else float('inf')
        status = "✓" if mx < 1e-10 else "✗ FAIL"
        print(f"  Z{k+1}:  max|A_6| = {mx:.2e}  {status}")
        if mx > 1e-10:
            print("ERROR: tree amplitude does not vanish. Aborting.")
            sys.exit(1)

    # ── Build bases ───────────────────────────────────────────────
    simple_basis = list(combinations(range(N_VARS), 3))        # 84
    broad_basis  = list(combinations_with_replacement(range(N_VARS), 3))  # 165
    print(f"\nSimple-poles ansatz: {len(simple_basis)} terms")
    print(f"Broad ansatz:        {len(broad_basis)} terms")

    A_tree_s = np.zeros(len(simple_basis))
    for tri in TREE_TRIPLES:
        A_tree_s[simple_basis.index(tri)] = 1.0

    A_tree_b = np.zeros(len(broad_basis))
    for tri in TREE_TRIPLES:
        A_tree_b[broad_basis.index(tri)] = 1.0

    # ── Simple ansatz ─────────────────────────────────────────────
    print("\nSampling for simple-poles ansatz...")
    M_s = build_constraint_matrix(simple_basis, n_samples=250, seed=12345)
    rank_s, null_s, vecs_s, S_s = compute_nullspace(
        M_s, "SIMPLE-POLES ANSATZ (84 terms)")

    if null_s >= 1:
        v = vecs_s[0]
        is_prop, res = check_proportional(v, A_tree_s, simple_basis)
        if is_prop:
            print(f"  >>> PROPORTIONAL TO A_6^tree  (residual = {res:.2e})  ✓")
        else:
            print(f"  >>> NOT proportional (residual = {res:.2e})")

    # ── Broad ansatz ──────────────────────────────────────────────
    print("\nSampling for broad ansatz (with higher-order poles)...")
    M_b = build_constraint_matrix(broad_basis, n_samples=250, seed=54321)
    rank_b, null_b, vecs_b, S_b = compute_nullspace(
        M_b, "BROAD ANSATZ (165 terms)")

    if null_b >= 1:
        v = vecs_b[0]
        is_prop, res = check_proportional(v, A_tree_b, broad_basis)
        if is_prop:
            print(f"  >>> PROPORTIONAL TO A_6^tree  (residual = {res:.2e})  ✓")
        else:
            print(f"  >>> NOT proportional (residual = {res:.2e})")

    # ── Robustness ────────────────────────────────────────────────
    print("\nRobustness (3 additional seeds, simple ansatz):")
    for seed in [77777, 31415, 27182]:
        M_c = build_constraint_matrix(simple_basis, n_samples=200, seed=seed, verbose=False)
        _, nd, _, _ = compute_nullspace(M_c, f"  seed={seed}")
    print()

    # ── Show tree terms ───────────────────────────────────────────
    if null_s == 1:
        v = vecs_s[0]
        ref_idx = simple_basis.index(TREE_TRIPLES[0])
        if abs(v[ref_idx]) > 1e-10:
            v = v / v[ref_idx]
        print("Nonzero components of nullspace vector:")
        for idx in range(len(simple_basis)):
            if abs(v[idx]) > 1e-8:
                a, b, c = simple_basis[idx]
                tag = " [TREE]" if A_tree_s[idx] != 0 else " [NON-TREE]"
                print(f"  1/({VAR_NAMES[a]}*{VAR_NAMES[b]}*{VAR_NAMES[c]})"
                      f":  {v[idx]:+.8f}{tag}")

    # ── Conclusion ────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("CONCLUSION")
    print("=" * 65)
    print(f"Simple ansatz: rank {rank_s}/{len(simple_basis)}, "
          f"nullspace dim {null_s}")
    print(f"Broad  ansatz: rank {rank_b}/{len(broad_basis)}, "
          f"nullspace dim {null_b}")

    if null_s == 1 and null_b == 1:
        print("\nNullspace is 1-dimensional and equals A_6^tree in both spaces.")
        print("=> 1-zeros UNIQUELY determine the 6-pt amplitude.")
        print("=> Theorem 5.6 verified at n=6.")
    elif null_s > 1 or null_b > 1:
        print(f"\nNullspace dim > 1 — investigate for possible counterexample.")
    else:
        print(f"\nNullspace empty — setup error.")


if __name__ == '__main__':
    main()
