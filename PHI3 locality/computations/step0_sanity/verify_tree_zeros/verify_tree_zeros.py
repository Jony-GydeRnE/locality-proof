#!/usr/bin/env python3
"""
verify_tree_zeros.py — Symbolically verify that the 5-pt phi^3 tree amplitude
vanishes on each of the 5 cyclic 1-zero loci.

PURPOSE
-------
This is a direct, self-contained verification of the "hidden zeros" property:
the phi^3 tree amplitude A_5 = 0 on each Z_k for k = 1,...,5.

This corresponds to the forward direction of Theorem 1 in the paper. The harder
converse direction (A_5 is the UNIQUE such function) is established by the nullspace
computation in full_ansatz_n5.py.

KINEMATICS
----------
Five planar Mandelstam invariants: X13, X14, X24, X25, X35
(These correspond to X_{ij} = (p_i + ... + p_{j-1})^2 in the paper's notation.)

The phi^3 tree amplitude for 5 external legs is the sum over the 5 non-crossing
triangulations of a convex pentagon with vertices 1,2,3,4,5:

  A_5 = 1/(X13*X14) + 1/(X13*X35) + 1/(X14*X24) + 1/(X24*X25) + 1/(X25*X35)

1-ZERO LOCI
-----------
Each locus Z_k imposes 2 linear constraints on the 5 variables, leaving a
3-dimensional subspace. The constraints come from the general formula (Eq. 1.2
in the paper):

  Z1: X13 = -X25,      X14 = X24 - X25
  Z2: X13 = -X25+X35,  X24 = X25 - X35
  Z3: X13 = X14+X35,   X24 = -X35
  Z4: X14 = -X35,      X24 = X25 - X35
  Z5: X13 = -X25+X35,  X14 = -X25

HOW TO RUN
----------
  python3 verify_tree_zeros.py

Requires: sympy (pip install sympy)
Runtime: < 1 second.
"""

from sympy import symbols, simplify, factor, expand, together

# ── Variables ─────────────────────────────────────────────────────────────────
# Using X_ij notation to match the paper's planar Mandelstam conventions.
X13, X14, X24, X25, X35 = symbols('X13 X14 X24 X25 X35')

# ── The 5-pt phi^3 tree amplitude ─────────────────────────────────────────────
# Sum over all 5 non-crossing triangulations of the convex pentagon (1,2,3,4,5).
# Each term 1/(X_ij * X_kl) corresponds to a unique triangulation.
# See Eq. (1.3) / Eq. (3) of the paper.
A5 = (1/(X13*X14)    # triangulation using diagonals 1-3 and 1-4
    + 1/(X13*X35)    # triangulation using diagonals 1-3 and 3-5
    + 1/(X14*X24)    # triangulation using diagonals 1-4 and 2-4
    + 1/(X24*X25)    # triangulation using diagonals 2-4 and 2-5
    + 1/(X25*X35))   # triangulation using diagonals 2-5 and 3-5

print("=" * 60)
print("Verification: A_5^tree vanishes on all 5 cyclic 1-zero loci")
print("=" * 60)
print()
print("The phi^3 tree amplitude:")
print("  A_5 = 1/(X13*X14) + 1/(X13*X35) + 1/(X14*X24)")
print("       + 1/(X24*X25) + 1/(X25*X35)")
print()
print("The 5 cyclic 1-zero loci (free variables are the unlisted ones):")
print("  Z1: X13 = -X25,       X14 = X24 - X25     (free: X24, X25, X35)")
print("  Z2: X13 = -X25 + X35, X24 = X25 - X35     (free: X14, X25, X35)")
print("  Z3: X13 = X14 + X35,  X24 = -X35          (free: X14, X25, X35)")
print("  Z4: X14 = -X35,       X24 = X25 - X35     (free: X13, X25, X35)")
print("  Z5: X13 = -X25 + X35, X14 = -X25          (free: X24, X25, X35)")
print()

# ── Define the 5 loci as substitution dictionaries ────────────────────────────
# Each locus solves two of the five variables in terms of the remaining three.
loci = {
    1: {X13: -X25,          X14: X24 - X25},
    2: {X13: -X25 + X35,    X24: X25 - X35},
    3: {X13: X14 + X35,     X24: -X35},
    4: {X14: -X35,          X24: X25 - X35},
    5: {X13: -X25 + X35,    X14: -X25},
}

# ── Verify A_5 = 0 on each locus ──────────────────────────────────────────────
print("Verification results:")
print("-" * 60)

all_pass = True
for k in range(1, 6):
    sub = loci[k]

    # Substitute the locus constraints into A_5
    A5_on_locus = A5.subs(sub)

    # Simplify to a canonical form (combine into single fraction, reduce)
    A5_simplified = simplify(A5_on_locus)

    is_zero = (A5_simplified == 0)
    status = "PASS ✓" if is_zero else "FAIL ✗"
    all_pass = all_pass and is_zero

    # Show the substitution and result
    sub_str = ",  ".join(f"{str(k)} → {str(v)}" for k, v in sub.items())
    print(f"\nLocus Z{k}:")
    print(f"  Substitute: {sub_str}")
    print(f"  A_5 |_{{Z{k}}} = {A5_simplified}")
    print(f"  Status: {status}")

    # If nonzero, show the factored form to help debug
    if not is_zero:
        factored = factor(A5_on_locus)
        print(f"  Factored form: {factored}")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)
if all_pass:
    print("All 5 loci verified: A_5^tree = 0 on Z_k for k = 1,...,5.")
    print()
    print("This confirms the 'hidden zeros' property of the phi^3 tree amplitude.")
    print("The forward direction of Theorem 1 is verified symbolically.")
    print()
    print("The converse (A_5 is the UNIQUE function satisfying this) is")
    print("established by the nullspace computation in full_ansatz_n5.py.")
else:
    print("ERROR: Some loci failed! Check the locus definitions.")
