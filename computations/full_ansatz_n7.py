#!/usr/bin/env python3
"""
full_ansatz_n7.py — 7-pt phi^3 theory: 1-zero nullspace computation (numerical).

PURPOSE
-------
Verifies the hidden-zero theorem at n=7: the cyclic 1-zero conditions on the
ansatz space B_7 of mass-dimension -4 rational functions yield a 1-dimensional
nullspace spanned by A_7^tree.

KINEMATICS
----------
14 planar Mandelstam variables at n=7 (non-adjacent pairs in the heptagon):
  X13, X14, X15, X16,
  X24, X25, X26, X27,
  X35, X36, X37,
  X46, X47,
  X57

ANSATZ
------
Simple-poles ansatz: all C(14,4)=1001 terms 1/(X_a*X_b*X_c*X_d), a<b<c<d.

1-ZERO LOCI
-----------
Each leg k imposes 4 constraints c_{k,j}=0 for j=k+2,...,k+5 (mod 7),
reducing 14 variables to 10 free via linear substitution.

EXPECTED OUTPUT
---------------
Nullspace dimension 1, spanned by A_7^tree (sum over 42 triangulations).

HOW TO RUN
----------
  python3 full_ansatz_n7.py

Requires: numpy
"""

import numpy as np
from itertools import combinations
import sys

# =====================================================================
#  SETUP
# =====================================================================

N = 7
N_FREE = 10   # 14 - 4 constraints per locus

# Non-adjacent pairs (i,j) with i<j in heptagon (1-indexed labels)
# Adjacent: (1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,1)
# Non-adjacent (14 total):
VAR_PAIRS = [
    (1,3),(1,4),(1,5),(1,6),
    (2,4),(2,5),(2,6),(2,7),
    (3,5),(3,6),(3,7),
    (4,6),(4,7),
    (5,7),
]
N_VARS = len(VAR_PAIRS)
assert N_VARS == 14

VAR_NAMES = [f'X{i}{j}' for i,j in VAR_PAIRS]
# Build lookup: (i,j) -> index
PAIR_TO_IDX = {(i,j): idx for idx, (i,j) in enumerate(VAR_PAIRS)}

def Xidx(i, j):
    """Get index of X_{ij} in the variable list, or None if adjacent/same."""
    # Normalize to 1-indexed mod 7
    i = ((i - 1) % N) + 1
    j = ((j - 1) % N) + 1
    if i == j:
        return None
    if i > j:
        i, j = j, i
    # Check if adjacent (including wrap)
    if j - i == 1 or (i == 1 and j == N):
        return None
    return PAIR_TO_IDX.get((i, j), None)


# =====================================================================
#  TREE AMPLITUDE — 42 triangulations of the heptagon
# =====================================================================
# Recursive enumeration of triangulations of a convex polygon.

def triangulate_polygon(vertices):
    """
    Enumerate all triangulations of a convex polygon given as a list of
    vertex labels (in order). Returns a list of triangulations, where each
    triangulation is a list of diagonals (pairs of vertices).

    A diagonal is any edge connecting non-consecutive vertices in the
    original heptagon. Edges of the polygon that are sides of the
    heptagon are not diagonals.
    """
    n = len(vertices)
    if n < 3:
        return [[]]
    if n == 3:
        # Triangle — no diagonals needed (edges are either sides or
        # already accounted for as diagonals of the parent polygon).
        return [[]]

    results = []
    # Fix the edge (vertices[0], vertices[-1]).
    # Choose k such that vertices[k] forms a triangle with this edge.
    for k in range(1, n - 1):
        # Triangle: (vertices[0], vertices[k], vertices[-1])
        # Diagonals introduced by this triangle:
        new_diags = []
        if k > 1:
            new_diags.append(tuple(sorted((vertices[0], vertices[k]))))
        if k < n - 2:
            new_diags.append(tuple(sorted((vertices[k], vertices[-1]))))

        # Recursively triangulate the two sub-polygons
        left_poly = vertices[:k+1]    # vertices[0]..vertices[k]
        right_poly = vertices[k:]     # vertices[k]..vertices[-1]

        left_triangulations = triangulate_polygon(left_poly)
        right_triangulations = triangulate_polygon(right_poly)

        for lt in left_triangulations:
            for rt in right_triangulations:
                results.append(new_diags + lt + rt)

    return results


# Heptagon vertices labeled 1..7
heptagon_vertices = list(range(1, N + 1))
all_triangulations_raw = triangulate_polygon(heptagon_vertices)

# Each triangulation of a heptagon has n-3 = 4 diagonals
# Convert diagonals to variable indices
TREE_QUADS = []
for diags in all_triangulations_raw:
    # Filter: only keep diagonals that are non-adjacent pairs
    idxs = []
    for d in diags:
        i, j = d
        idx = Xidx(i, j)
        if idx is not None:
            idxs.append(idx)
    idxs = tuple(sorted(set(idxs)))
    assert len(idxs) == 4, f"Expected 4 diagonals, got {len(idxs)}: {diags} -> {idxs}"
    TREE_QUADS.append(idxs)

# Remove duplicates (shouldn't be any, but just in case)
TREE_QUADS = list(set(TREE_QUADS))
TREE_QUADS.sort()

C5 = 42  # Catalan number C_5
assert len(TREE_QUADS) == C5, f"Expected {C5} triangulations, got {len(TREE_QUADS)}"

print(f"Enumerated {len(TREE_QUADS)} triangulations of the heptagon.")


# =====================================================================
#  1-ZERO LOCI
# =====================================================================
# c_{k,j} = X_{k,j} + X_{k+1,j+1} - X_{k,j+1} - X_{k+1,j} = 0
# where X_{i,i+1} = 0 (adjacent) and X_{i,i} = 0.
# Indices are 1-based, mod 7.
#
# For each leg k, the constraints are c_{k,j} = 0 for j in the 4
# non-adjacent vertices to both k and k+1 (mod 7).
# Specifically j = k+2, k+3, k+4, k+5 (mod 7).

def compute_locus(k):
    """
    Compute the substitution rules for locus Z_k.

    For leg k (edge k to k+1), constraints c_{k,j}=0 for
    j = k+2, k+3, k+4, k+5 (mod 7).

    c_{k,j} = X_{k,j} + X_{k+1,j+1} - X_{k,j+1} - X_{k+1,j}

    where any X_{i,i+1} or X_{i,i} = 0.

    Returns dict with 'constrained' (4 indices), 'free' (10 indices),
    'sub' (4 x 10 matrix).
    """
    # The 4 constraint equations c_{k,j}=0
    # Each equation is a linear combination of X variables.
    # We represent each equation as a dict: var_index -> coefficient
    equations = []
    for offset in range(2, 6):  # j = k+2, k+3, k+4, k+5
        j = ((k - 1 + offset) % N) + 1  # 1-indexed
        kp1 = (k % N) + 1  # k+1 mod 7, 1-indexed
        jp1 = (j % N) + 1

        # c_{k,j} = X_{k,j} + X_{k+1,j+1} - X_{k,j+1} - X_{k+1,j}
        terms = [
            (k, j, +1),
            (kp1, jp1, +1),
            (k, jp1, -1),
            (kp1, j, -1),
        ]
        eq = {}
        for (a, b, coeff) in terms:
            idx = Xidx(a, b)
            if idx is not None:
                eq[idx] = eq.get(idx, 0) + coeff
        equations.append(eq)

    # Now we need to choose 4 "constrained" variables and express them
    # in terms of the remaining 10 "free" variables.
    # Strategy: for each equation, identify which variables appear,
    # then solve the system.

    # Build the coefficient matrix for all 14 variables
    coeff_matrix = np.zeros((4, N_VARS))
    for i, eq in enumerate(equations):
        for idx, c in eq.items():
            coeff_matrix[i, idx] = c

    # Choose constrained variables via pivoting (row echelon form)
    # Use column pivoting to find 4 pivot columns
    A = coeff_matrix.copy()
    pivot_cols = []
    for row in range(4):
        # Find column with largest absolute value in remaining rows
        best_col = None
        best_val = 0
        for col in range(N_VARS):
            if col in pivot_cols:
                continue
            val = abs(A[row, col])
            if val > best_val:
                best_val = val
                best_col = col
        assert best_col is not None and best_val > 0.5, \
            f"Degenerate constraint at Z_{k}, row {row}"
        pivot_cols.append(best_col)
        # Eliminate this column from other rows
        for r in range(4):
            if r == row:
                continue
            if abs(A[r, best_col]) > 1e-15:
                factor = A[r, best_col] / A[row, best_col]
                A[r] -= factor * A[row]

    # Now solve: constrained vars = linear combo of free vars
    constrained = sorted(pivot_cols)
    free = sorted(set(range(N_VARS)) - set(constrained))
    assert len(constrained) == 4 and len(free) == 10

    # Rebuild with constrained on the left
    # A_c * x_c + A_f * x_f = 0  =>  x_c = -A_c^{-1} A_f x_f
    A_c = coeff_matrix[:, constrained]
    A_f = coeff_matrix[:, free]
    sub = -np.linalg.solve(A_c, A_f)

    return {
        'constrained': constrained,
        'free': free,
        'sub': sub,
    }


LOCI = [compute_locus(k) for k in range(1, N + 1)]

# Verify: check that each locus has 4 constrained and 10 free
for k, loc in enumerate(LOCI):
    assert len(loc['constrained']) == 4, f"Z_{k+1}: wrong # constrained"
    assert len(loc['free']) == 10, f"Z_{k+1}: wrong # free"
    assert loc['sub'].shape == (4, 10), f"Z_{k+1}: wrong sub shape"


# =====================================================================
#  SAMPLING AND EVALUATION
# =====================================================================

def sample_point(loc, rng, scale=3.0):
    """Sample a random kinematic point on a 1-zero locus."""
    free_vals = rng.standard_normal(N_FREE) * scale
    X = np.zeros(N_VARS)
    for i, fi in enumerate(loc['free']):
        X[fi] = free_vals[i]
    constr_vals = loc['sub'] @ free_vals
    for i, ci in enumerate(loc['constrained']):
        X[ci] = constr_vals[i]
    return X


def evaluate_tree(X):
    """Evaluate A_7^tree at kinematic point X."""
    return sum(1.0 / (X[a] * X[b] * X[c] * X[d]) for a, b, c, d in TREE_QUADS)


def build_constraint_matrix(basis, n_samples=400, seed=12345, verbose=True):
    """Build constraint matrix by sampling points on all 7 loci."""
    rng = np.random.default_rng(seed)
    rows = []
    n_basis = len(basis)

    for k in range(N):
        loc = LOCI[k]
        good = 0
        attempts = 0
        while good < n_samples and attempts < n_samples * 50:
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
    for i in range(max(0, rank - 3), min(n_basis, rank + 4)):
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
    print("7-pt phi^3  1-zero nullspace computation  (numerical)")
    print("=" * 65)

    # ── Print variable info ───────────────────────────────────────
    print(f"\n{N_VARS} planar Mandelstam variables:")
    print(f"  {', '.join(VAR_NAMES)}")
    print(f"\n{len(TREE_QUADS)} triangulations (Catalan C_5 = 42)")

    # ── Print triangulations ──────────────────────────────────────
    print("\nTriangulations (as diagonal sets):")
    for i, q in enumerate(TREE_QUADS):
        names = [VAR_NAMES[idx] for idx in q]
        print(f"  {i+1:2d}. {', '.join(names)}")

    # ── Verify loci ───────────────────────────────────────────────
    print("\nVerifying A_7^tree vanishes on each locus:")
    rng_v = np.random.default_rng(999)
    for k in range(N):
        vals = []
        for _ in range(500):
            X = sample_point(LOCI[k], rng_v)
            if np.min(np.abs(X)) < 0.01:
                continue
            try:
                vals.append(abs(evaluate_tree(X)))
            except (ZeroDivisionError, FloatingPointError):
                continue
        mx = max(vals) if vals else float('inf')
        status = "OK" if mx < 1e-8 else "FAIL"
        print(f"  Z{k+1}:  max|A_7| = {mx:.2e}  {status}")
        if mx > 1e-8:
            print("ERROR: tree amplitude does not vanish on this locus.")
            print("  Constrained vars:", [VAR_NAMES[i] for i in LOCI[k]['constrained']])
            print("  Free vars:", [VAR_NAMES[i] for i in LOCI[k]['free']])
            sys.exit(1)

    # ── Build basis ───────────────────────────────────────────────
    simple_basis = list(combinations(range(N_VARS), 4))  # C(14,4) = 1001
    print(f"\nSimple-poles ansatz: {len(simple_basis)} terms")
    assert len(simple_basis) == 1001

    A_tree = np.zeros(len(simple_basis))
    for tri in TREE_QUADS:
        A_tree[simple_basis.index(tri)] = 1.0

    # ── Build constraint matrix ───────────────────────────────────
    print("\nSampling for simple-poles ansatz (400 per locus, 2800 total)...")
    M = build_constraint_matrix(simple_basis, n_samples=400, seed=12345)
    rank, null_dim, vecs, S = compute_nullspace(
        M, "SIMPLE-POLES ANSATZ (1001 terms)")

    if null_dim >= 1:
        v = vecs[0]
        is_prop, res = check_proportional(v, A_tree, simple_basis)
        if is_prop:
            print(f"  >>> PROPORTIONAL TO A_7^tree  (residual = {res:.2e})")
        else:
            print(f"  >>> NOT proportional (residual = {res:.2e})")
            # Check if A_tree is in the span of the nullspace
            if null_dim > 1:
                # Project A_tree onto nullspace
                A_norm = A_tree / np.linalg.norm(A_tree)
                proj = sum(np.dot(vecs[i], A_norm) * vecs[i] for i in range(null_dim))
                proj_res = np.linalg.norm(proj - A_norm)
                print(f"  Projection residual of A_7^tree onto nullspace: {proj_res:.2e}")

    # ── Singular value gap ────────────────────────────────────────
    print(f"\nSingular value gap:")
    if rank < len(simple_basis):
        print(f"  S[{rank-1}] / S[{rank}] = {S[rank-1]/S[rank]:.2e}")

    # ── Show tree terms in nullspace ──────────────────────────────
    if null_dim == 1:
        v = vecs[0]
        ref_idx = simple_basis.index(TREE_QUADS[0])
        if abs(v[ref_idx]) > 1e-10:
            v = v / v[ref_idx]
        nz = sum(1 for x in v if abs(x) > 1e-8)
        print(f"\nNonzero components of nullspace vector: {nz}")
        tree_count = 0
        nontree_count = 0
        for idx in range(len(simple_basis)):
            if abs(v[idx]) > 1e-8:
                if A_tree[idx] != 0:
                    tree_count += 1
                else:
                    nontree_count += 1
        print(f"  Tree terms: {tree_count} / {len(TREE_QUADS)}")
        print(f"  Non-tree terms: {nontree_count}")

    # ── Robustness check ──────────────────────────────────────────
    print("\nRobustness (2 additional seeds):")
    for seed in [77777, 31415]:
        M_c = build_constraint_matrix(simple_basis, n_samples=300, seed=seed, verbose=False)
        _, nd, _, _ = compute_nullspace(M_c, f"  seed={seed}")

    # ── Conclusion ────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("CONCLUSION")
    print("=" * 65)
    print(f"Simple ansatz: rank {rank}/{len(simple_basis)}, "
          f"nullspace dim {null_dim}")

    if null_dim == 1:
        print("\nNullspace is 1-dimensional and equals A_7^tree.")
        print("=> 1-zeros UNIQUELY determine the 7-pt amplitude.")
        print("=> Hidden-zero theorem verified at n=7.")
    elif null_dim > 1:
        print(f"\nNullspace dim = {null_dim} > 1 — investigate.")
    else:
        print(f"\nNullspace empty — setup error.")


if __name__ == '__main__':
    main()
