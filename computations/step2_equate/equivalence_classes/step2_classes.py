"""
================================================================================
STEP-2 EQUIVALENCE-CLASS ANALYSIS
"After Step 1 kills the obvious non-locals, what equality chains does Step 2
produce among the survivors?  And do all triangulations land in ONE class?"
================================================================================

What Step 2 is (in plain words)
-------------------------------
At every cyclic zone Z_{r, r+2} the 1-zero conditions force the bare relation

    X_{r+1, r-1}  =  - X_{r, r+2}               (e.g., X_{25} = -X_{13} at Z_{1,3})

In Step 1 we sent the special X_{r,r+2} to infinity; this killed coefficients.
In Step 2 we keep X_{r,r+2} FINITE and read what the bare relation forces on
the coefficients themselves.

Naively: the monomial 1/X_{bare} on the zone equals -1/X_{special}, so a
multiset M containing the bare contributes the same rational form (up to a
sign) as the multiset M' obtained by swapping the bare for the special.
Demanding B|_zone = 0 then equates those coefficients: a_M = a_{M'}.

By doing this at every zone (cyclic shifts), we get a CHAIN of equalities
that partitions the ansatz coefficients into EQUIVALENCE CLASSES.  Within
each class all coefficients are equal.  Combined with Step 1, three things
can happen to a class:

  (i)  PURE SURVIVING:  every member of the class is uncaught by Step 1.
        Then all members have one common (potentially nonzero) value.
  (ii) PURE KILLED:  every member is killed by Step 1.  Trivial (= 0 = 0).
  (iii) MIXED:  some members survive, some are killed.  The equality forces
        a_M = 0 for ALL members.  In effect Step 2 KILLS the survivors in
        this class via their connection to a Step-1-killed sibling.

Unitarity statement:  if all C_{n-2} triangulations end up in ONE pure-
surviving class, then all triangulation coefficients are equal -> the tree
amplitude is uniquely determined up to one overall scale.

A flag worth knowing about
--------------------------
The CLEAN equality a_M = a_{M'} only follows from the bare relation when no
"messy" terms in B|_zone contribute to the same rational form (e.g.,
contributions from substituted chords' partial fractions).  The full Step-2
system is a multi-term linear system over all coefficients.

**Two modes** are computed and reported below, to bracket the truth:

  - PERMISSIVE mode:  every bare<->special swap (bare in M, special not in
    M, no substitutes from this zone in M) gives an edge.  This rule
    over-equates because it ignores the multi-term partial-fraction
    contributions from substituted chords in OTHER multisets.  At n=5
    (where unitarity is known to hold from full_ansatz_n5.py), permissive
    mode incorrectly groups every triangulation into a class with
    Step-1-killed siblings -- a clear signal of overcounting.

  - CONSERVATIVE mode (the meaningful one):  only count swap edges whose
    BOTH endpoints are Step-1 survivors.  This is a LOWER BOUND on the
    true equivalence (we never spuriously chain to a killed multiset), so
    if the conservative rule already puts all triangulations in ONE class,
    that's strong evidence for unitarity that doesn't depend on the
    multi-term machinery.

Usage
-----
    python3 computations/step2_classes.py            # prompts for n
    python3 computations/step2_classes.py 7
    python3 computations/step2_classes.py 8

Recommended cap: n <= 8 in pure Python (n=9 with ~900k multisets and an
equivalence-graph build takes a few minutes; cap for sanity).
"""

import argparse
import os
import sys
import time
from itertools import combinations_with_replacement


# ============================================================================
# 1.  COMBINATORICS OF CHORDS AND ZONES (same as step1_uncaught.py / step1_dual.py)
# ============================================================================

def normalize(i, j, n):
    """Canonical name (a, b) of the chord between vertices i and j, with
    a < b in {1,...,n}.  Indices may be given mod n."""
    i = ((i - 1) % n) + 1
    j = ((j - 1) % n) + 1
    if i > j:
        i, j = j, i
    return (i, j)


def all_chords(n):
    """Every chord of the n-gon: pairs (i, j), 1 <= i < j <= n, that are
    NOT polygon edges.  Total = n(n-3)/2."""
    return [(i, j) for i in range(1, n + 1)
                   for j in range(i + 1, n + 1)
                   if j - i >= 2 and not (i == 1 and j == n)]


def zone_structure(r, n):
    """Special, bare, and (companion, substitute) pairs at zone Z_{r, r+2}.
    See step1_uncaught.py for details."""
    special = normalize(r,     r + 2, n)
    bare    = normalize(r + 1, r - 1, n)
    pairs = []
    for offset in range(3, n - 1):
        k = ((r - 1 + offset) % n) + 1
        pairs.append((normalize(r, k, n), normalize(r + 1, k, n)))
    return special, bare, pairs


def killable_at_zone(ms_set, structure):
    """Step-1 kill criterion at this zone.  Returns True iff the multiset
    can be killed by sending the special to infinity at this zone."""
    special, bare, pairs = structure
    if special in ms_set or bare in ms_set:
        return False
    for c1, c2 in pairs:
        if c1 in ms_set and c2 in ms_set:
            return False
    return True


def chords_cross(c1, c2):
    a, b = c1; c, d = c2
    if {a, b} & {c, d}:
        return False
    inside = lambda v: a < v < b
    return inside(c) != inside(d)


def is_triangulation(ms):
    """A size-(n-3) multiset is a triangulation iff distinct + no crossings."""
    if len(set(ms)) != len(ms):
        return False
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            if chords_cross(ms[i], ms[j]):
                return False
    return True


# ============================================================================
# 2.  STEP-2 BARE<->SPECIAL EQUIVALENCE
# ============================================================================

def bare_swap_neighbors(M, n, structures):
    """
    Yield every multiset M' obtained from M by a SINGLE bare<->special swap
    at some zone, under the strict condition that the swap produces a
    clean equality a_M = a_{M'}.

    Conditions for the clean equality at zone Z_{r, r+2}:
      - M contains the bare chord X_{r+1, r-1}
      - M does NOT contain the special chord X_{r, r+2}
        (avoids "swap creates two specials" multi-term cases)
      - M contains no substituted chord X_{r+1, k}
        (otherwise that chord's partial-fraction expansion contributes to
        the same rational form as our swap, making the equation multi-term)

    F^triangle chords and companion chords X_{r, k} are OK to be present.

    The reverse direction (M containing special, swap to bare) is handled
    automatically because we iterate over ALL multisets and union-find is
    symmetric.  We need only do bare->special here.
    """
    M_list = list(M)
    M_set  = set(M)
    for r in range(1, n + 1):
        special, bare, pairs = structures[r - 1]
        if bare    not in M_set:    continue
        if special     in M_set:    continue
        substitutes = {sub for _, sub in pairs}
        if any(sub in M_set for sub in substitutes): continue
        # apply swap: replace one occurrence of bare with special
        i = M_list.index(bare)
        M_new = tuple(sorted(M_list[:i] + [special] + M_list[i + 1:]))
        yield M_new


def equiv_classes(items, n, restrict_to=None):
    """Compute equivalence classes via union-find on `items`.

    `restrict_to`, if given, is a set; an edge M~M' is added only if BOTH
    M and M' are in `restrict_to`. This implements CONSERVATIVE mode
    (only equate Step-1 survivors).  When `restrict_to` is None we use
    PERMISSIVE mode (any swap, even to a Step-1-killed multiset).
    """
    structures = [zone_structure(r, n) for r in range(1, n + 1)]
    parent = {M: M for M in items}

    def find(M):
        root = M
        while parent[root] != root:
            root = parent[root]
        while parent[M] != root:
            parent[M], M = root, parent[M]
        return root

    def union(M1, M2):
        r1, r2 = find(M1), find(M2)
        if r1 != r2:
            parent[r1] = r2

    items_set = set(items)
    for M in items:
        for M_neighbor in bare_swap_neighbors(M, n, structures):
            if M_neighbor not in items_set:
                continue
            if restrict_to is not None:
                if M not in restrict_to or M_neighbor not in restrict_to:
                    continue
            union(M, M_neighbor)

    classes = {}
    for M in items:
        rep = find(M)
        classes.setdefault(rep, []).append(M)
    return list(classes.values())


# ============================================================================
# 3.  CYCLIC SHIFT ANALYSIS
#     (helps decide whether classes are "shifts of triangulations")
# ============================================================================

def cyclic_shift(M, n, k=1):
    """Apply the Z_n action: every chord (i,j) -> ((i+k) mod n, (j+k) mod n)."""
    return tuple(sorted(normalize(i + k, j + k, n) for (i, j) in M))


def cyclic_orbit(M, n):
    """Orbit of M under the cyclic group Z_n (set of multisets reached by
    rotation)."""
    return {cyclic_shift(M, n, k) for k in range(n)}


def cyclic_orbit_reps(items, n):
    """(Used for individual multisets, not classes.)  Pick one representative
    per cyclic-shift orbit from a set of multisets."""
    items_set = set(items)
    seen = set()
    reps, sizes = [], []
    for M in items:
        if M in seen:
            continue
        orbit = cyclic_orbit(M, n) & items_set
        seen.update(orbit)
        reps.append(M)
        sizes.append(len(orbit))
    return reps, sizes


def class_orbit_decomposition(classes, n):
    """Group equivalence classes into orbits under the cyclic Z_n action.
    Two classes C and C' are in the same orbit iff some cyclic shift maps
    C (as a set of multisets) onto C'."""
    class_sets = [frozenset(cls) for cls in classes]
    parent = list(range(len(class_sets)))

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    # Build a lookup so we can match shifted classes quickly
    set_to_index = {s: i for i, s in enumerate(class_sets)}
    for i, C in enumerate(class_sets):
        for k in range(1, n):
            shifted = frozenset(cyclic_shift(M, n, k) for M in C)
            j = set_to_index.get(shifted)
            if j is not None:
                union(i, j)

    orbits = {}
    for i in range(len(class_sets)):
        rep = find(i)
        orbits.setdefault(rep, []).append(i)
    return list(orbits.values())


# ============================================================================
# 4.  MAIN
# ============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Step-2 bare<->special equivalence-class analysis."
    )
    parser.add_argument("n", nargs="?", type=int, default=None,
                        help="Number of legs (>=4, recommended <=8).")
    args = parser.parse_args()
    if args.n is None:
        args.n = int(input("Enter n (>= 4, recommended <= 8): ").strip())
    return args.n


def main():
    n = parse_args()
    if n < 4:
        sys.exit(f"n = {n} too small; need n >= 4.")
    N = n - 3

    chords = all_chords(n)
    structures = [zone_structure(r, n) for r in range(1, n + 1)]
    all_ms = list(combinations_with_replacement(chords, N))
    print(f"\nn = {n}, multiset size N = {N}, |chords| = {len(chords)}")
    print(f"Total size-{N} multisets: {len(all_ms)}\n")

    # -------- Step 1 survivors ----------
    print("Computing Step-1 survivors (all n zones)...")
    t0 = time.time()
    survivors = set()
    for M in all_ms:
        ms_set = set(M)
        if not any(killable_at_zone(ms_set, s) for s in structures):
            survivors.add(M)
    print(f"  ({time.time()-t0:.2f} s)")
    tris    = [m for m in survivors if is_triangulation(m)]
    nontris = [m for m in survivors if not is_triangulation(m)]
    print(f"Step-1 survivors: {len(survivors):5d}  =  "
          f"{len(tris):4d} triangulations  +  "
          f"{len(nontris):4d} non-triangulations\n")

    # ============== PERMISSIVE MODE (full multiset space) ==============
    print("=" * 60)
    print("PERMISSIVE mode (over-equates; for diagnostic only)")
    print("=" * 60)
    t0 = time.time()
    permissive_classes = equiv_classes(all_ms, n, restrict_to=None)
    print(f"  ({time.time()-t0:.2f} s)")

    surv_perm   = [c for c in permissive_classes if all(m in survivors for m in c)]
    kill_perm   = [c for c in permissive_classes if all(m not in survivors for m in c)]
    mixed_perm  = [c for c in permissive_classes
                   if any(m in survivors for m in c)
                   and any(m not in survivors for m in c)]
    print(f"Total classes: {len(permissive_classes)}")
    print(f"  pure-surviving: {len(surv_perm)}   "
          f"pure-killed: {len(kill_perm)}   MIXED: {len(mixed_perm)}")
    if mixed_perm:
        n_tris_in_mixed = sum(1 for c in mixed_perm
                              for m in c if m in tris and m in survivors)
        print(f"  -> {n_tris_in_mixed} triangulation(s) appear in MIXED classes,")
        print(f"     which (under permissive rule) would falsely force them to 0.")
        print(f"     The permissive rule is OVER-equating; see CONSERVATIVE mode below.")
    print()

    # ============== CONSERVATIVE MODE (survivors only) ==============
    print("=" * 60)
    print("CONSERVATIVE mode (survivors only; meaningful answer)")
    print("=" * 60)
    t0 = time.time()
    survivors_list = sorted(survivors)
    cons_classes = equiv_classes(survivors_list, n, restrict_to=survivors)
    print(f"  ({time.time()-t0:.2f} s)")
    print(f"Step-2 equivalence classes among Step-1 survivors: {len(cons_classes)}")
    print()

    # size distribution
    from collections import Counter
    size_counts = Counter(len(cls) for cls in cons_classes)
    print("Class size distribution:")
    for sz in sorted(size_counts):
        print(f"  size {sz:3d}: {size_counts[sz]} class(es)")
    print()

    # triangulation distribution
    tri_set = set(tris)
    classes_containing_tri = [cls for cls in cons_classes
                              if any(m in tri_set for m in cls)]
    print(f"Triangulations land in {len(classes_containing_tri)} class(es).")
    if len(classes_containing_tri) == 1:
        cls = classes_containing_tri[0]
        tris_in = [m for m in cls if m in tri_set]
        non_in  = [m for m in cls if m not in tri_set]
        print(f"  *** UNITARITY (conservative Step-2): "
              f"all {len(tris_in)} triangulations in ONE class of size {len(cls)}.")
        if non_in:
            print(f"      That class also contains {len(non_in)} non-triangulation survivor(s),")
            print(f"      which the conservative rule equates to the triangulation value.")
    else:
        print(f"  Triangulations split across multiple classes:")
        for i, cls in enumerate(classes_containing_tri):
            tris_in    = [m for m in cls if m in tri_set]
            non_tri_in = [m for m in cls if m not in tri_set]
            print(f"   class {i+1}: size {len(cls)}, "
                  f"{len(tris_in)} tri + {len(non_tri_in)} non-tri")
    print()

    # non-triangulation survivor placement
    if nontris:
        nontri_in_tri_class = sum(1 for m in nontris
                                  if any(m in cls for cls in classes_containing_tri))
        nontri_in_own = len(nontris) - nontri_in_tri_class
        print(f"Non-triangulation Step-1 survivors ({len(nontris)} total):")
        print(f"  - In a triangulation-class (forced equal to triangulation value): "
              f"{nontri_in_tri_class}")
        print(f"  - In their own (singleton or no-triangulation) classes:           "
              f"{nontri_in_own}")
        print()

    # cyclic-orbit decomposition of the conservative classes themselves
    orbits = class_orbit_decomposition(cons_classes, n)
    orbit_size_counts = Counter(len(orb) for orb in orbits)
    print("Cyclic-shift orbits of the conservative classes:")
    print(f"  {len(cons_classes)} classes  ->  {len(orbits)} cyclic orbit(s).")
    for sz in sorted(orbit_size_counts):
        print(f"    {orbit_size_counts[sz]} orbit(s) of size {sz}")
    print()

    # interpretation summary
    print("=" * 60)
    print("READING THE OUTPUT")
    print("=" * 60)
    print("- Permissive shows what naive bare<->special swaps give if you")
    print("  trust every swap; it OVER-equates because it ignores multi-term")
    print("  contributions from substituted-chord partial fractions, and so")
    print("  spuriously links Step-1 survivors to Step-1-killed siblings.")
    print("- Conservative restricts to swaps where BOTH endpoints already")
    print("  survive Step 1; it is a LOWER bound on the true Step-2")
    print("  equivalence (any further indirect chain through a non-survivor")
    print("  intermediate is suppressed).")
    print("- The full Step-2 system (with multi-term equations) sits between")
    print("  these two bounds.  Existing full_ansatz_n5/6/7.py confirm the")
    print("  full system has nullspace dim = 1 -> all triangulation")
    print("  coefficients ARE in fact equal -> unitarity holds.")
    print("- So if conservative shows triangulations split into >1 class,")
    print("  the missing equivalences come from indirect chains through")
    print("  non-survivor intermediates that the conservative rule blocks.")
    print()


if __name__ == "__main__":
    main()
