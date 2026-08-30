#!/usr/bin/env python3
"""
build_full_system.py
====================

Empirically computes the rank of the *full* linear system on the multiset
coefficients a_M coming from the vanishing of B_n on every cyclic hidden
k-zero Z^(k)_r, ranging over all k in {1,..,n-3} and all cyclic rotations
r in {1,..,n}.

Why this script exists
----------------------
Paper Part II classifies the *kills* on each Z^(k) (the rate-vector
A/B/C criteria). Paper Part III classifies the *Step-2 coefficient
relations* on each Z^(k). Both are projections of one bigger object:

    B_n |_{Z^(k)_r}  =  0   as a Laurent polynomial in the free coordinates

When we equate every Laurent monomial coefficient of that vanishing to
zero, we get a linear functional on the a_M. The union of all such
functionals over all (k, r) is the *full constraint system*.

The all-n locality + unitarity conjecture says: this system has
nullspace dimension 1, spanned by the "triangulation indicator"
(a_M = 1 if M is a triangulation, 0 otherwise) -- equivalently every
non-triangulation coefficient vanishes and all triangulation
coefficients equal a single common scalar.

This script tests the conjecture at n in {5, 6, 7} by sampling random
points on each Z^(k)_r, building constraint rows, and computing the
rank modulo a large prime via incremental Gaussian elimination.

Method
------
For each (k, r):
  1. Parametrize Z^(k)_r by random rational values of the free coords
     (special X_{1,k+2}, offsets X_{i,k+2} for i=2..k, column-1 chords
     X_{1,l} for l=k+3..n-1, and the born-free chords). All canonical
     labels are cyclically shifted by r.
  2. Compute every chord X_{i,j} on Z^(k)_r via Part II Lemma 1.
  3. Each random point gives one constraint:
        sum_M (1 / prod_{c in M} X_c(point)) * a_M = 0
  4. Convert to mod-PRIME and feed into an incremental row-reduction.

We keep sampling until the rank stabilises (no growth across a full
sweep over all (k, r) pairs) or hits T(n) - 1, where T(n) is the
number of multisets.

Outputs land in ./outputs/n{N}_results.{md, json}.
"""
import os
import json
import random
import time
from fractions import Fraction
from itertools import combinations_with_replacement, combinations
import numpy as np

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

# NTT-friendly prime: 998244353 = 119 * 2^23 + 1. Fits comfortably in int64
# arithmetic (since 998244353^2 < 2^63).
PRIME = 998244353

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
os.makedirs(OUT_DIR, exist_ok=True)


# ----------------------------------------------------------------------
# Chord and multiset enumeration
# ----------------------------------------------------------------------

def all_chords(n):
    """All chord variables (i,j) with 1<=i<j<=n, not a leg.

    Legs are X_{i,i+1} and X_{1,n} (the two sides of the polygon at each
    vertex). Chord count = n(n-3)/2."""
    return [(i, j) for i in range(1, n + 1) for j in range(i + 2, n + 1)
            if not (i == 1 and j == n)]


def all_multisets(n):
    """All multisets of size n-3 over chords; count = C(N+n-4, n-3)."""
    return list(combinations_with_replacement(all_chords(n), n - 3))


def chords_cross(c1, c2):
    """Two chords (a,b) and (c,d) cross iff a<c<b<d or c<a<d<b."""
    a, b = sorted(c1)
    c, d = sorted(c2)
    return (a < c < b < d) or (c < a < d < b)


def is_triangulation(M, n):
    """A multiset M is a triangulation iff it is a SET (no repetitions)
    of n-3 chords with no pairwise crossing."""
    if len(set(M)) != n - 3:
        return False
    for i, j in combinations(range(n - 3), 2):
        if chords_cross(M[i], M[j]):
            return False
    return True


# ----------------------------------------------------------------------
# Cyclic shift helpers
# ----------------------------------------------------------------------

def shift_vertex(v, r, n):
    """Cyclically shift vertex v by r-1 mod n (1-indexed). r=1 is the identity."""
    return ((v - 1 + r - 1) % n) + 1


def shift_chord(c, r, n):
    """Cyclic shift of chord (a,b); returns sorted pair."""
    a = shift_vertex(c[0], r, n)
    b = shift_vertex(c[1], r, n)
    if a > b:
        a, b = b, a
    return (a, b)


# ----------------------------------------------------------------------
# Random point on the k-zero constraint surface Z^(k)_r
# ----------------------------------------------------------------------

def sample_point(n, k, r, rng):
    """Sample a random rational point on Z^(k)_r.

    Returns: dict {(a,b): Fraction value} populated for every chord (a,b).
    Returns None if the sampled point is singular (a chord came out zero),
    so the caller should regenerate.

    The k-zero based at vertices {r, r+1, ..., r+k} (cyclic mod n,
    1-indexed) has free coords (after the cyclic shift by r):

      - special:     X_{sigma(1), sigma(k+2)}
      - offsets:     X_{sigma(i), sigma(k+2)} for i = 2..k
      - column-1:    X_{sigma(1), sigma(l)} for l = k+3..n-1
      - born-free:   chords whose CANONICAL labels (a,b) have both
                     endpoints in {1..k+1} or both in {k+2..n}.

    The remaining chords are determined by Part II Lemma 1:
        X_{sigma(k+1), sigma(l)} = X_{sigma(1), sigma(l)} - special
        X_{sigma(k+1), sigma(n)} = - special                       (bare)
        X_{sigma(i),  sigma(j)}  = X_{sigma(k+1), sigma(j)} + X_{sigma(i), sigma(k+2)}
    """
    def randval():
        # Random nonzero rational. Small range keeps numerators tiny so the
        # mod-p reductions stay sane.
        while True:
            v = rng.randint(-25, 25)
            if v != 0:
                return Fraction(v)

    def s(c):
        return shift_chord(c, r, n)

    X = {}

    # ---- Free coordinates ----
    chord_special = s((1, k + 2))
    special_val = randval()
    X[chord_special] = special_val

    for i in range(2, k + 1):
        cc = (min(i, k + 2), max(i, k + 2))
        X[s(cc)] = randval()

    for l in range(k + 3, n):
        X[s((1, l))] = randval()

    # Born-free chords: those internal to {1..k+1} or to {k+2..n} (canonical).
    for c in all_chords(n):
        a, b = c
        if (1 <= a < b <= k + 1) or (k + 2 <= a < b <= n):
            cr = s(c)
            if cr not in X:
                X[cr] = randval()

    # ---- Determined chords (Part II Lemma 1, in rotated form) ----
    # Backbone: X_{sigma(k+1), sigma(l)} = X_{sigma(1), sigma(l)} - special
    for l in range(k + 3, n):
        c_kp1_l = s((min(k + 1, l), max(k + 1, l)))
        c_1_l   = s((1, l))
        X[c_kp1_l] = X[c_1_l] - special_val

    # Bare: X_{sigma(k+1), sigma(n)} = - special
    c_bare = s((min(k + 1, n), max(k + 1, n)))
    if c_bare not in X:
        X[c_bare] = -special_val

    # Interior bridges: X_{sigma(i), sigma(j)} = X_{sigma(k+1), sigma(j)} + X_{sigma(i), sigma(k+2)}
    for i in range(2, k + 1):
        c_i_kp2 = s((min(i, k + 2), max(i, k + 2)))
        t_i = X[c_i_kp2]
        for j in range(k + 3, n + 1):
            c_kp1_j = s((min(k + 1, j), max(k + 1, j)))
            c_i_j   = s((min(i, j), max(i, j)))
            if c_i_j not in X:
                X[c_i_j] = X[c_kp1_j] + t_i

    # Sanity / singularity check
    chords = all_chords(n)
    for c in chords:
        if c not in X:
            return None  # should not happen; bug in substitution
        if X[c] == 0:
            return None  # singular: regenerate
    return X


# ----------------------------------------------------------------------
# Constraint row construction (mod PRIME)
# ----------------------------------------------------------------------

def build_row(X, multisets, p):
    """Given chord values X (dict of Fractions) and a list of multisets,
    return the constraint row mod p:

        row[idx_M] = (prod_{c in M} X[c])^{-1} mod p

    Each random point gives the equation sum_M row[idx_M] * a_M = 0.

    Returns None if any chord value mod p is zero (then this row is
    pathological; the caller should skip)."""
    row = np.zeros(len(multisets), dtype=np.int64)
    # Precompute each chord's mod-p value once
    chord_mod = {}
    for c, val in X.items():
        num = int(val.numerator) % p
        den = int(val.denominator) % p
        if num == 0 or den == 0:
            return None
        chord_mod[c] = (num * pow(den, p - 2, p)) % p
        if chord_mod[c] == 0:
            return None
    for idx, M in enumerate(multisets):
        prod = 1
        for c in M:
            prod = (prod * chord_mod[c]) % p
        row[idx] = pow(prod, p - 2, p)  # 1 / prod  mod p
    return row


# ----------------------------------------------------------------------
# Incremental rank over GF(p)
# ----------------------------------------------------------------------

class ModRank:
    """Incremental rank computation over GF(p) using row reduction.

    Maintains a list of pivot rows in reduced echelon form. Adding a new
    row reduces it against existing pivots; if nonzero, it becomes a new
    pivot. Rank = number of pivots.

    Note: this is O(rank * T) per row, fine for T up to a few thousand.
    For larger T we'd want a smarter data structure (sparse, block
    Wiedemann, ...). For our n <= 7 scope, this is plenty fast."""

    def __init__(self, T, p):
        self.T = T
        self.p = p
        self.pivots = []  # list of (pivot_col, row_array)

    def add_row(self, row):
        """Reduce `row` against existing pivots; add it as new pivot if
        nonzero. Returns True iff rank increased."""
        p = self.p
        row = row.copy() % p
        for pc, prow in self.pivots:
            if row[pc] != 0:
                scale = int(row[pc])
                row = (row - scale * prow) % p
        nz = np.flatnonzero(row)
        if len(nz) == 0:
            return False
        pc_new = int(nz[0])
        inv = pow(int(row[pc_new]), p - 2, p)
        row = (row * inv) % p
        # Eliminate this new pivot column from existing pivot rows (RREF)
        for i, (pc, prow) in enumerate(self.pivots):
            if prow[pc_new] != 0:
                scale = int(prow[pc_new])
                self.pivots[i] = (pc, (prow - scale * row) % p)
        self.pivots.append((pc_new, row))
        return True

    @property
    def rank(self):
        return len(self.pivots)


# ----------------------------------------------------------------------
# Nullspace verification
# ----------------------------------------------------------------------

def vector_in_nullspace(incr, v):
    """Is the vector v (length T) in the null space, i.e. each pivot row
    has dot product zero mod p with v?"""
    p = incr.p
    for pc, prow in incr.pivots:
        if int((prow * v).sum() % p) != 0:
            return False
    return True


def nullspace_basis(incr, T):
    """Compute a basis for the null space from the RREF pivot rows.

    For each free (non-pivot) column c, build a basis vector with
    v[c]=1 and v[pivot_col]=-row[c] for each pivot row.

    Returns a list of numpy int64 arrays of length T, all in the
    null space."""
    p = incr.p
    pivot_cols = set(pc for pc, _ in incr.pivots)
    free_cols = [c for c in range(T) if c not in pivot_cols]
    basis = []
    for fc in free_cols:
        v = np.zeros(T, dtype=np.int64)
        v[fc] = 1
        for pc, prow in incr.pivots:
            v[pc] = (-int(prow[fc])) % p
        basis.append(v)
    return basis


# ----------------------------------------------------------------------
# Main pipeline per n
# ----------------------------------------------------------------------

def run_for_n(n, seed=42, max_samples=None, verbose=True):
    """Sample random points across all (k, r), build the constraint
    matrix mod PRIME, and compute the rank.

    Returns (incr, T, samples_used, elapsed_seconds, kr_pairs)."""
    multisets = all_multisets(n)
    T = len(multisets)
    if max_samples is None:
        max_samples = max(8000, T * 4)
    target_rank = T - 1  # what we expect if the conjecture holds

    rng = random.Random(seed)
    incr = ModRank(T, PRIME)

    # (k, r) pairs: k = 1..n-3, r = 1..n
    kr_pairs = [(k, r) for k in range(1, n - 2) for r in range(1, n + 1)]
    if verbose:
        print(f"n = {n}: T(n) = {T}, |kr| = {len(kr_pairs)}, "
              f"target_rank = {target_rank}")

    t0 = time.time()
    samples_used = 0
    no_growth_passes = 0
    prev_pass_rank = -1

    while samples_used < max_samples:
        for (k, r) in kr_pairs:
            X = sample_point(n, k, r, rng)
            if X is None:
                continue
            row = build_row(X, multisets, PRIME)
            if row is None:
                continue
            incr.add_row(row)
            samples_used += 1
            if incr.rank >= target_rank:
                break
        elapsed = time.time() - t0
        if verbose:
            print(f"  pass: samples={samples_used}, rank={incr.rank}, "
                  f"elapsed={elapsed:.1f}s")
        if incr.rank >= target_rank:
            break
        if incr.rank == prev_pass_rank:
            no_growth_passes += 1
            if no_growth_passes >= 3:
                if verbose:
                    print("  no rank growth in 3 passes; stopping.")
                break
        else:
            no_growth_passes = 0
        prev_pass_rank = incr.rank

    elapsed = time.time() - t0
    return incr, T, samples_used, elapsed, kr_pairs


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------

def make_report(n, incr, T, multisets, samples_used, elapsed, kr_pairs):
    """Build a markdown report + JSON summary of the run."""
    # Triangulation indicator vector
    triang_vec = np.zeros(T, dtype=np.int64)
    triang_indices = []
    for idx, M in enumerate(multisets):
        if is_triangulation(M, n):
            triang_vec[idx] = 1
            triang_indices.append(idx)
    triang_count = len(triang_indices)
    triang_in_null = vector_in_nullspace(incr, triang_vec)

    null_dim = T - incr.rank

    # Nullspace basis (only computed if dim is small, else skip for speed)
    null_basis_info = None
    if null_dim <= 5 and null_dim >= 1:
        basis = nullspace_basis(incr, T)
        null_basis_info = []
        for v in basis:
            # Find indices where v is nonzero, plus distinguish triang vs non-triang
            nz = np.flatnonzero(v)
            nz_triang = [int(i) for i in nz if i in set(triang_indices)]
            nz_nontriang = [int(i) for i in nz if i not in set(triang_indices)]
            null_basis_info.append({
                "n_nonzero": int(len(nz)),
                "n_nonzero_triang": len(nz_triang),
                "n_nonzero_nontriang": len(nz_nontriang),
                "matches_triang_indicator": (
                    len(nz_nontriang) == 0
                    and len(nz_triang) == triang_count
                    and all(v[i] == v[triang_indices[0]] for i in triang_indices)
                ),
            })

    summary = {
        "n": n,
        "T": T,
        "C_{n-2}": triang_count,
        "kr_pairs": len(kr_pairs),
        "samples_used": samples_used,
        "rank": incr.rank,
        "nullspace_dim": null_dim,
        "triangulation_indicator_in_nullspace": bool(triang_in_null),
        "elapsed_seconds": round(elapsed, 2),
        "prime": PRIME,
        "nullspace_basis": null_basis_info,
    }

    # Markdown
    lines = []
    lines.append(f"# Full constraint-system rank at n = {n}\n")
    lines.append(f"- T(n) (multisets of size {n-3}): **{T}**")
    lines.append(f"- C_{{n-2}} (triangulations): **{triang_count}**")
    lines.append(f"- (k, r) pairs sampled: **{len(kr_pairs)}** "
                 f"(k = 1..{n-3}, r = 1..{n})")
    lines.append(f"- Constraints sampled: **{samples_used}**")
    lines.append(f"- **Rank** of constraint matrix: **{incr.rank}**")
    lines.append(f"- **Nullspace dimension**: **{null_dim}**")
    lines.append(f"- Triangulation-indicator vector in nullspace: "
                 f"**{triang_in_null}**")
    lines.append(f"- Prime: {PRIME}")
    lines.append(f"- Elapsed: {elapsed:.1f} s")
    lines.append("")

    if null_dim == 1 and triang_in_null:
        lines.append("## Verdict")
        lines.append("")
        lines.append("**Cokernel is 1-dim, spanned by the triangulation indicator.**")
        lines.append("")
        lines.append("Every non-triangulation coefficient a_M is forced to zero, and")
        lines.append(f"every triangulation coefficient equals one common scalar.")
        lines.append("This is locality + unitarity at this n from the FULL k-zero")
        lines.append("system -- strictly stronger than what Part I proves (k=1 only).")
    elif null_dim == 1 and not triang_in_null:
        lines.append("## Verdict")
        lines.append("")
        lines.append("**Nullspace is 1-dim, but NOT spanned by triangulation indicator.**")
        lines.append("This would be a striking obstruction -- locality fails as stated.")
    elif null_dim > 1:
        lines.append("## Verdict")
        lines.append("")
        lines.append(f"**Nullspace is {null_dim}-dim** (expected 1).")
        if triang_in_null:
            lines.append("The triangulation indicator IS in the nullspace, but there")
            lines.append("are other kernel directions: the full k-zero system at this n")
            lines.append("does not close to 1-dim by sampling alone. Either more sampling")
            lines.append("is needed (rank may not have converged), or the conjecture")
            lines.append("genuinely requires additional non-rate-vector input.")
        else:
            lines.append("The triangulation indicator is NOT in the nullspace.")
        if null_basis_info:
            lines.append("")
            lines.append("### Nullspace basis vectors")
            for idx, info in enumerate(null_basis_info):
                lines.append(
                    f"- vector {idx}: {info['n_nonzero']} nonzero "
                    f"({info['n_nonzero_triang']} triang + "
                    f"{info['n_nonzero_nontriang']} non-triang); "
                    f"matches triang indicator: {info['matches_triang_indicator']}"
                )
    elif null_dim == 0:
        lines.append("## Verdict")
        lines.append("")
        lines.append("**Rank = T(n); nullspace is trivial.** Over-determined relative")
        lines.append("to expectation -- this is likely a sampling artifact (random")
        lines.append("rationals giving a spurious extra constraint).")

    return "\n".join(lines), summary


# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

def main(ns=(5, 6, 7), seed=42):
    """Run for each n in `ns`; write outputs/n{n}_results.{md,json}."""
    overall = {}
    for n in ns:
        print(f"\n=========== n = {n} ===========")
        incr, T, samples, elapsed, kr_pairs = run_for_n(n, seed=seed)
        multisets = all_multisets(n)
        report_md, summary = make_report(
            n, incr, T, multisets, samples, elapsed, kr_pairs
        )

        md_path = os.path.join(OUT_DIR, f"n{n}_results.md")
        with open(md_path, "w") as f:
            f.write(report_md + "\n")
        json_path = os.path.join(OUT_DIR, f"n{n}_results.json")
        with open(json_path, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"  wrote {md_path}")
        print(f"  wrote {json_path}")
        overall[n] = summary

    # Combined summary
    combined_path = os.path.join(OUT_DIR, "summary.json")
    with open(combined_path, "w") as f:
        json.dump(overall, f, indent=2)
    print(f"\nwrote {combined_path}")
    return overall


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        ns = tuple(int(x) for x in sys.argv[1:])
    else:
        ns = (5, 6, 7)
    main(ns)
