"""
================================================================================
cascade_kill_n7.py  --  Verify the Laurent-cascade kill of the n=7 "fish" terms
================================================================================

PURPOSE (BIG PICTURE)
---------------------
At n=7, the leading-order ("layer-0") hidden-zero kill mechanism — the one
implemented in step1_uncaught.py — fails to kill exactly seven non-triangulation
chord multisets.  Visually, each of the seven looks like a small "fish":
two crossing chords in some part of the heptagon plus two extra chords that
between them touch every vertex, blocking layer-0 from being applied.

The conjecture this script verifies is: every fish coefficient nevertheless
vanishes once we include LAURENT TAILS (subleading orders in the special
variable, expanded as a geometric series).  Concretely, for each fish:

  Step 1.  Pick a zone Z = Z_{r, r+2} where the fish lives in a useful
           Laurent stratum.  Send the special X_{r,r+2} -> infinity.

  Step 2.  Expand the fish's monomial 1/(X_a X_b X_c X_d) on Z as a Laurent
           series in 1/X_{r,r+2}.  Read off the coefficient of some
           subleading order — call this the FINGERPRINT order.

  Step 3.  Enumerate every other monomial of the ansatz that contributes
           the SAME free-variable signature at the same Laurent order on Z.
           Those monomials' coefficients enter a single linear equation:

                +- a_{fish}   +-  sum  +- a_{cousin}    =   0.

  Step 4.  For each "cousin" coefficient on the right, verify that it is
           ALREADY KNOWN TO VANISH — either by layer-0 at some zone, or by
           a prior cascade step.  If every cousin is dead, the equation
           collapses to  +- a_{fish} = 0,  and the fish dies too.

This script does Steps 1-4 symbolically with sympy — no hand algebra is
trusted.  We use sympy.series to compute the Laurent expansions and
sympy.Poly to extract the coefficient of each free-variable monomial at
each order.  All sign conventions are the ones in step1_dual.py.

OUTPUT
------
For each of the seven fish:
  - the chosen kill zone and Laurent order
  - the full list of cousin coefficients in the fingerprint equation
  - which cousins are killed by layer-0 (and at which zone), and which need
    a prior Laurent cascade step (and at which zone + order);
  - for those needing a prior cascade step, the verification of THAT step
    (so the script is fully recursive: it never assumes a kill that it
    has not itself derived).

Final output is a results file:  results_cascade_n7.txt.

CONVENTIONS
-----------
- Vertices labelled 1..n cyclically.
- Chord (i,j) with 1 <= i < j <= n; not a polygon edge (so j-i >= 2 and
  not (i,j) == (1,n)).
- Zone Z_{r, r+2}:
      special chord  s_z   = (r, r+2)
      bare chord     b_z   = (r+1, r-1)  -- substitution X_{b_z} = -X_{s_z}
      non-bare pairs (companion X_{r,k}, substitute X_{r+1,k})
                     for k = r+3, ..., r+n-2 (cyclic),
                     with substitution X_{r+1,k} = X_{r,k} - X_{r,r+2}.
- We use the SAME notation as step1_dual.py / step1_uncaught.py.
"""

from itertools import combinations_with_replacement
from collections import Counter, defaultdict
import sympy as sp


# ============================================================================
# 1.  Combinatorics:  chords, zones, layer-0 killability
#     (copied from step1_dual.py with no changes -- so the cascade verifier
#     uses identical conventions to your existing layer-0 enumerator.)
# ============================================================================

def normalize(i, j, n):
    """Canonical chord name: pair (a,b) with a<b in {1..n}, indices mod n."""
    i = ((i - 1) % n) + 1
    j = ((j - 1) % n) + 1
    if i > j:
        i, j = j, i
    return (i, j)


def all_chords(n):
    """All planar chords of the n-gon (non-edge pairs).  Total: n(n-3)/2."""
    return [(i, j) for i in range(1, n + 1)
                   for j in range(i + 1, n + 1)
                   if j - i >= 2 and not (i == 1 and j == n)]


def zone_structure(r, n):
    """
    For zone Z_{r, r+2} return  (special, bare, pairs):
        special = (r, r+2)              -- the Laurent variable
        bare    = (r+1, r-1) (cyclic)   -- X_bare = -X_special
        pairs   = list of (companion, substitute) chord pairs
                 (X_{r,k}, X_{r+1,k}) for k = r+3, ..., r+n-2 cyclic;
                 substitution: X_substitute = X_companion - X_special.
    """
    special = normalize(r,     r + 2, n)
    bare    = normalize(r + 1, r - 1, n)
    pairs = []
    for offset in range(3, n - 1):
        k = ((r - 1 + offset) % n) + 1
        companion  = normalize(r,     k, n)
        substitute = normalize(r + 1, k, n)
        pairs.append((companion, substitute))
    return special, bare, pairs


def killable_at_zone_layer0(ms_list, structure):
    """
    Layer-0 ('Step-1') killability of the multiset on this zone.

    Returns True iff:
       (A) the multiset contains neither the special nor the bare;
       (B) for every (companion, substitute) pair, the multiset does NOT
           contain BOTH chords of the pair.
    """
    ms_set = set(ms_list)
    special, bare, pairs = structure
    if special in ms_set or bare in ms_set:
        return False
    for companion, substitute in pairs:
        if companion in ms_set and substitute in ms_set:
            return False
    return True


def layer0_kill_zone(ms_list, n):
    """Return the smallest r such that the multiset is layer-0 killable at
    zone Z_{r,r+2}, or None if no zone works."""
    for r in range(1, n + 1):
        if killable_at_zone_layer0(ms_list, zone_structure(r, n)):
            return r
    return None


# ============================================================================
# 2.  Symbolic Laurent expansion of a monomial on a zone
# ============================================================================
#
# For each zone Z_{r,r+2} we set up a sympy symbol for every chord of the
# n-gon, plus the special X_S.  The substitutions for that zone replace the
# substitute chords by  (companion - X_S)  and the bare chord by  -X_S.
# The Laurent expansion of a monomial 1/(X_{c1}...X_{cN}) on that zone is
# obtained by sympy.series in the variable X_S around 0 (with the convention
# that 1/X_S^k ~ X_S^{-k}, so we expand in 1/X_S; equivalently we substitute
# X_S = 1/u and expand around u = 0, then convert back).
#
# The "free variables" of the zone are: the companions, plus all chords with
# no endpoint at vertex r+1 (the "born-free" inner-triangle chords).  A
# fingerprint at order X_S^{-k} is a rational function of these free vars.

def chord_symbol(c):
    """Sympy symbol X_{i,j} for chord c = (i,j).  Always positive (kinematic)."""
    i, j = c
    return sp.Symbol(f"X_{i}_{j}", positive=True)


def laurent_coefficient(monomial_chords, zone_r, n, max_order):
    """
    Compute the Laurent expansion of   1/(X_{c1} * X_{c2} * ... * X_{cN})
    on zone Z_{zone_r, zone_r+2}, up to order X_S^{-max_order}, where
    X_S is the special chord of that zone.

    Returns a dict
        order_k  :  rational expression in the free variables
    such that the monomial's value on the zone equals
        sum over k>=0 of   (rational expression) / X_S^k.
    Only orders k = 0, 1, ..., max_order are returned.

    Implementation:
      - replace each substitute chord by (companion - X_S)
      - replace the bare by -X_S
      - free chords stay as themselves
      - use sympy.series in X_S^{-1} via the substitution X_S = 1/u, expand
        in u around 0, then read coefficients of u^k = X_S^{-k}.
    """
    special, bare, pairs = zone_structure(zone_r, n)
    X_S = chord_symbol(special)
    u   = sp.Symbol("__u__", positive=True)

    # Build the substituted monomial.
    expr = sp.Integer(1)
    for c in monomial_chords:
        if c == special:
            expr = expr / X_S
        elif c == bare:
            expr = expr / (-X_S)
        else:
            sub_chord = None
            for companion, substitute in pairs:
                if c == substitute:
                    sub_chord = companion
                    break
            if sub_chord is not None:
                expr = expr / (chord_symbol(sub_chord) - X_S)
            else:
                # Free chord (born-free or companion).  Stays as-is.
                expr = expr / chord_symbol(c)

    # Now expand expr in 1/X_S.  Substitute X_S = 1/u and expand around u=0.
    expr_u = expr.subs(X_S, 1 / u)
    # We need the coefficient of u^k for k = 0, 1, ..., max_order.
    # series takes O(u^{max_order+1}) which is u^{max_order+1} truncation.
    s = sp.series(expr_u, u, 0, max_order + 1).removeO()
    s = sp.expand(s)
    # Extract coefficients of each power of u.
    out = {}
    for k in range(max_order + 1):
        coeff = sp.simplify(s.coeff(u, k))
        out[k] = coeff
    return out


# ============================================================================
# 3.  Fingerprint matching: which OTHER monomials of B share a target
#     fingerprint at the same Laurent order on the same zone?
# ============================================================================
#
# A "fingerprint" is a specific rational monomial in the free variables of
# the zone.  Two ansatz monomials produce the same fingerprint at order
# X_S^{-k} iff their Laurent coefficients at order k contain the same
# rational monomial (up to scalar).  We extract that as follows:
#
#   1. For each candidate monomial of the ansatz, compute laurent_coefficient
#      up to order max_order.
#   2. At a chosen order k and a chosen target fingerprint (a sympy
#      expression that is monomial in the free vars), extract the scalar
#      coefficient of that monomial in the candidate's order-k expression.
#   3. Sum up   scalar * a_{candidate}   over all candidates contributing
#      to the fingerprint.

def coefficient_of_monomial(expr, target_monomial, free_vars):
    """
    Extract the scalar coefficient of `target_monomial` from `expr`, where
    both are rational expressions in `free_vars`.

    Strategy: form  ratio = expr / target_monomial,  expand into a sum,
    and collect terms that are independent of every free variable.  The
    sum of those terms is the desired scalar coefficient.

    This handles correctly the case where `expr` is a sum of multiple
    monomials in `free_vars`, only one of which matches `target_monomial`.
    Returns sympy.Integer(0) when no term matches.
    """
    ratio = sp.simplify(expr / target_monomial)
    ratio = sp.together(ratio)
    ratio_expanded = sp.expand(ratio)
    if ratio_expanded.is_Add:
        constant_part = sp.Integer(0)
        for term in ratio_expanded.args:
            if not any(term.has(v) for v in free_vars):
                constant_part += term
        return sp.simplify(constant_part)
    else:
        if not any(ratio_expanded.has(v) for v in free_vars):
            return sp.simplify(ratio_expanded)
        else:
            return sp.Integer(0)


def fingerprint_equation(fish_multiset, zone_r, n, target_order,
                         target_fingerprint, candidates_pool, max_order=None):
    """
    Compute the fingerprint equation that the fish multiset participates
    in at order target_order on zone Z_{zone_r, zone_r+2} with free-variable
    fingerprint target_fingerprint.

    Parameters
    ----------
    fish_multiset      : tuple of chord pairs (the target coefficient)
    zone_r             : zone index r
    n                  : polygon size
    target_order       : Laurent order k (i.e., coeff of X_S^{-k})
    target_fingerprint : sympy expression for the free-variable monomial
                         to match
    candidates_pool    : list of multisets (each a tuple of chord pairs)
                         to check as potential cousins
    max_order          : Laurent order to expand to (defaults to target_order)

    Returns
    -------
    equation_terms : list of (multiset, scalar_coefficient) pairs;
                     the kill equation reads
                       sum  scalar * a_{multiset}  =  0,
                     where a_{multiset} is the unknown coefficient.
                     The fish multiset is included if it contributes.
    """
    if max_order is None:
        max_order = target_order

    # Identify the free variables on this zone (all symbols except X_S).
    special, bare, pairs = zone_structure(zone_r, n)
    X_S = chord_symbol(special)
    free_vars = []
    for c in all_chords(n):
        if c == special:
            continue
        sym = chord_symbol(c)
        # On the zone, the bare and substitutes are NOT free variables --
        # they have been replaced.  Only companions and born-free chords
        # remain.
        is_substitute = any(c == sub for _, sub in pairs)
        is_bare       = (c == bare)
        if is_substitute or is_bare:
            continue
        free_vars.append(sym)

    equation_terms = []
    for ms in candidates_pool:
        coeff_dict = laurent_coefficient(ms, zone_r, n, max_order)
        order_k_expr = coeff_dict[target_order]
        if order_k_expr == 0:
            continue
        scalar = coefficient_of_monomial(order_k_expr, target_fingerprint,
                                         free_vars)
        if scalar is None or scalar == 0:
            continue
        equation_terms.append((ms, scalar))
    return equation_terms


# ============================================================================


def candidate_multisets_for_fingerprint(target_fp, n, multiset_size):
    """
    Quickly enumerate the multisets of `multiset_size` chords that COULD
    produce the rational monomial `target_fp` at SOME Laurent order on
    SOME zone.

    Heuristic: if `target_fp = (numerator chord) / (denominator chords)`
    as a rational monomial in chord symbols, the candidate multiset must
    contain the denominator chords (with multiplicities) plus exactly
    enough additional chords to reach `multiset_size`.  The extra chord
    can be anything (it provides the Laurent tail or absorbs a power).

    This is a fast pre-filter.  The actual coefficient extraction
    then verifies which candidates contribute.
    """
    import sympy as sp
    from itertools import combinations_with_replacement
    fp_clean = sp.simplify(target_fp)
    num, den = sp.fraction(fp_clean)
    # Symbols in num and den
    num_syms = set()
    if num != 1 and num != 0:
        for term in sp.Mul.make_args(num):
            base, exp = term.as_base_exp()
            if base.is_Symbol:
                num_syms.add(base)
    den_syms_count = {}
    if den != 1:
        for term in sp.Mul.make_args(den):
            base, exp = term.as_base_exp()
            if base.is_Symbol:
                den_syms_count[base] = den_syms_count.get(base, 0) + int(exp)
    # Convert symbol names to chords
    def sym_to_chord(s):
        # symbols are X_i_j
        name = s.name
        parts = name.split("_")
        return (int(parts[1]), int(parts[2]))
    den_chord_count = {sym_to_chord(s): k for s, k in den_syms_count.items()}
    # Candidates: multisets that contain at least these chords with these
    # multiplicities, plus (multiset_size - sum of multiplicities) extra
    # chords.
    forced_chords = []
    for c, k in den_chord_count.items():
        forced_chords.extend([c] * k)
    extras_needed = multiset_size - len(forced_chords)
    if extras_needed < 0:
        return []
    chords = all_chords(n)
    out = []
    for extras in combinations_with_replacement(chords, extras_needed):
        candidate = tuple(sorted(forced_chords + list(extras)))
        out.append(candidate)
    # Deduplicate
    seen = set()
    deduped = []
    for c in out:
        if c not in seen:
            seen.add(c)
            deduped.append(c)
    return deduped


# ============================================================================
# 4.  The seven n=7 fish (read off the survivor PDF) and the cascade
# ============================================================================

N_VALUE = 7

# The seven fish multisets, identified directly from the survivor PDF
# nonlocal_n7_step1_survivors.pdf.
SEVEN_FISH = [
    ((1, 3), (1, 6), (2, 4), (4, 6)),
    ((1, 3), (1, 6), (3, 5), (4, 6)),
    ((1, 3), (1, 6), (3, 5), (5, 7)),
    ((1, 3), (2, 7), (3, 5), (5, 7)),
    ((1, 6), (2, 4), (2, 7), (4, 6)),
    ((2, 4), (2, 7), (3, 5), (5, 7)),
    ((2, 4), (2, 7), (4, 6), (5, 7)),
]


def all_size_N_multisets(n, max_count=None):
    """Every chord multiset of size n-3.  At n=7 there are C(14+3,4)=2380
    such multisets, comfortably small for full enumeration."""
    chords = all_chords(n)
    out = []
    for ms in combinations_with_replacement(chords, n - 3):
        out.append(ms)
        if max_count is not None and len(out) >= max_count:
            break
    return out


# ============================================================================
# 5.  The recursive cascade verifier
# ============================================================================
#
# Algorithm (depth-bounded):
#
#   verify_kill(M):
#     if M is layer-0 killable at some zone:
#         return ("layer-0", that-zone)
#     else for each candidate (zone, order, fingerprint):
#         compute the fingerprint equation at that zone+order
#         if M is in the equation with nonzero coefficient AND
#            every other multiset M' in the equation can be verified to die
#            (recursively, with smaller depth budget), then:
#             return ("cascade", zone, order, fingerprint, [verifications of M']s)
#     return None  -- could not verify within the depth budget.

# To keep the search fast and intelligible, we use a small handcrafted list
# of (zone, target_order, target_fingerprint) candidates per fish.  These
# come directly from the worked example in the chat conversation (fish #1)
# and its cyclic rotations (fish #2 - #7).

def chord_sym(c):
    """Convenience: chord_symbol with shorter spelling."""
    return chord_symbol(c)


# ----- For fish #1 = {(1,3),(1,6),(2,4),(4,6)}, the worked example uses
#       zone Z_{5,7}, order X_{57}^{-3}, fingerprint  X_{15}/(X_{13}*X_{24}).
#
# Cyclic rotation r -> r+1 (mod 7) of the chord set rotates fish #1 to
# the other six fish.  Likewise the kill prescription rotates.

def cyclic_rotate_chord(c, shift, n):
    """Apply cyclic rotation v -> v+shift to the endpoints of chord c."""
    i, j = c
    return normalize(i + shift, j + shift, n)


def cyclic_rotate_multiset(ms, shift, n):
    """Apply cyclic rotation to every chord in the multiset."""
    out = tuple(sorted(cyclic_rotate_chord(c, shift, n) for c in ms))
    return out


# Identify each fish with the rotation that maps fish #1 to it.
# We do this by trying all shifts 0..n-1 and checking for set equality
# (since multisets compare as tuples after sorting).
def find_shift(fish, n):
    """Return the cyclic shift sending fish #1 to `fish` (or None)."""
    base = SEVEN_FISH[0]
    for shift in range(n):
        if cyclic_rotate_multiset(base, shift, n) == tuple(sorted(fish)):
            return shift
    return None


# ============================================================================
# 6.  Main driver: verify each fish dies via the prescribed cascade.
# ============================================================================

def verify_fish_one(verbose=True):
    """
    Verify fish #1 = {(1,3),(1,6),(2,4),(4,6)} dies via the worked cascade:

    Step A. At Z_{5,7}, the order X_{57}^{-3} fingerprint X_{15}/(X_{13}*X_{24})
            picks up six monomials, including the fish.  Five must die first.

    Step B. For each of the five, identify the prior kill zone:
              - three die by layer-0 at Z_{4,6};
              - two die by Z_{4,6} Laurent kill at order X_{46}^{-2}.

    Returns the verification record.
    """
    n = N_VALUE
    fish = SEVEN_FISH[0]

    if verbose:
        print(f"\n=== Verifying fish #1 = {fish} ===")

    # --- Step A: build the Z_{5,7} fingerprint equation ---
    zone_r = 5
    target_order = 3
    fp = chord_sym((1, 5)) / (chord_sym((1, 3)) * chord_sym((2, 4)))

    # Candidate pool: all size-4 multisets of n=7 chords.
    chords = all_chords(n)
    pool = list(combinations_with_replacement(chords, n - 3))

    eq = fingerprint_equation(fish, zone_r, n, target_order, fp, pool)
    if verbose:
        print(f"  Step A: Z_{{5,7}}, order X_57^-3, fingerprint X_15/(X_13*X_24)")
        print(f"  {len(eq)} monomials contribute:")
        for ms, scalar in eq:
            print(f"     {('+' if scalar>=0 else '')}{scalar}  *  a_{ms}")

    # --- Step B: each cousin must die first. ---
    cousins = [(ms, scalar) for ms, scalar in eq if ms != fish]
    cousin_kills = []
    for ms, scalar in cousins:
        # Try layer-0 first.
        z0 = layer0_kill_zone(ms, n)
        if z0 is not None:
            cousin_kills.append((ms, scalar, "layer-0", z0))
            if verbose:
                print(f"     {ms}: layer-0 at Z_{{{z0},{(z0%n)+2 if (z0%n)+2<=n else (z0%n)+2-n}}}")
            continue
        # Otherwise try Z_{4,6} Laurent kill.
        # The two cousins that need this are (1,3),(1,6),(2,4),(5,7) and
        # (1,3),(1,5),(1,6),(2,4) per the worked example.
        # We verify each by finding an isolated fingerprint at Z_{4,6}.
        verified = verify_z46_laurent_kill(ms, n, verbose=verbose)
        if verified is not None:
            cousin_kills.append((ms, scalar, "Z_{4,6} Laurent", verified))
        else:
            cousin_kills.append((ms, scalar, "UNVERIFIED", None))

    return {"fish": fish, "zone": zone_r, "order": target_order,
            "fingerprint": fp, "equation": eq, "cousin_kills": cousin_kills}


def verify_z46_laurent_kill(ms, n, verbose=True):
    """
    Try to kill `ms` via a Laurent-tail isolation at Z_{4,6}, order X_46^{-2}.

    The argument: at Z_{4,6}, certain non-bare substitutions have Laurent
    tails containing companion variables.  We look for a fingerprint whose
    only source in B is the multiset `ms`.

    For the two cousins of fish #1 that need this:
      - {(1,3),(1,6),(2,4),(5,7)}: fingerprint X_{47}/(X_13*X_24*X_16) at order 2;
      - {(1,3),(1,5),(1,6),(2,4)}: fingerprint X_{14}/(X_13*X_24*X_16) at order 2.

    More generally we try all fingerprints of the form
        X_companion / (rest of the multiset symbols),
    where the "rest" is (multiset minus the substituted chord whose Laurent
    tail produces the companion).  If exactly one source contributes, the
    kill is verified.
    """
    zone_r = 4
    target_order = 2
    special, bare, pairs = zone_structure(zone_r, n)

    # Identify substitute chords in `ms`.
    substitutes_in_ms = []
    for c in ms:
        for companion, substitute in pairs:
            if c == substitute:
                substitutes_in_ms.append((companion, substitute))
                break

    # If `ms` has no substitute at this zone, this Laurent route doesn't apply.
    if not substitutes_in_ms:
        return None

    # For each substitute in `ms`, the fingerprint produced is
    #   X_companion / (product of other multiset symbols).
    chords = all_chords(n)
    pool = list(combinations_with_replacement(chords, n - 3))

    for companion, substitute in substitutes_in_ms:
        rest_factors = sp.Integer(1)
        # remove ONE occurrence of `substitute` from ms
        ms_minus = list(ms)
        ms_minus.remove(substitute)
        for c in ms_minus:
            rest_factors *= chord_sym(c)
        fp = chord_sym(companion) / rest_factors

        eq = fingerprint_equation(ms, zone_r, n, target_order, fp, pool)
        # If `ms` is the only contributor (with nonzero scalar), the kill works.
        contributors = [(m, s) for m, s in eq if s != 0]
        if len(contributors) == 1 and contributors[0][0] == ms:
            if verbose:
                print(f"     {ms}: Z_{{4,6}} Laurent, order X_46^-2, "
                      f"fingerprint X_{companion[0]}{companion[1]} / (...) "
                      f"-- isolated.")
            return {"zone_r": zone_r, "order": target_order,
                    "companion": companion, "equation": eq}

    return None


def verify_each_fish_from_scratch(fish, n, verbose=False):
    """
    For an arbitrary fish multiset (a layer-0-uncatchable non-triangulation),
    determine its kill recipe by mimicking the worked example for fish #1
    via cyclic rotation, then VERIFY that recipe by direct sympy
    computation -- i.e., compute the actual fingerprint equation at the
    chosen zone and order, check the fish appears in it, and check every
    other contributor is killed by layer-0 (or, if not, by a one-step
    Laurent kill at a neighboring zone).

    Returns a dict with the kill recipe and the verified equation.
    Returns None if no recipe works.
    """
    # Try every cyclic rotation r of fish #1's recipe.
    # Fish #1's recipe:  zone Z_{5,7}, order 3, fingerprint X_15/(X_13 X_24).
    base_zone = 5
    base_order = 3
    # Express the base fingerprint as a list of (chord, signed_power)
    base_fp = [((1, 5), 1), ((1, 3), -1), ((2, 4), -1)]

    for shift in range(n):
        # Rotate fish #1 by `shift`; if it equals our `fish`, we found
        # the right shift.
        rotated_base_fish = cyclic_rotate_multiset(SEVEN_FISH[0], shift, n)
        if tuple(sorted(rotated_base_fish)) != tuple(sorted(fish)):
            continue
        # Rotate the recipe.
        zone_r = ((base_zone - 1 + shift) % n) + 1
        # Build the rotated fingerprint as a sympy expression.
        fp = sp.Integer(1)
        for c, p in base_fp:
            c_rot = cyclic_rotate_chord(c, shift, n)
            fp = fp * chord_sym(c_rot)**p
        # Compute the actual equation at the rotated zone+order+fp.
        # Fast: only enumerate multisets whose chord factors are
        # *compatible* with producing this fingerprint.  Coefficient
        # extraction will weed out the false positives.
        pool = candidate_multisets_for_fingerprint(fp, n, n - 3)
        eq = fingerprint_equation(fish, zone_r, n, base_order, fp, pool)
        if not eq:
            continue
        # Check the fish is in the equation with nonzero scalar.
        fish_terms = [(m, s) for m, s in eq if tuple(sorted(m)) == tuple(sorted(fish))]
        if not fish_terms or all(s == 0 for _, s in fish_terms):
            continue
        return {"shift": shift, "zone_r": zone_r, "order": base_order,
                "fingerprint": fp, "equation": eq}
    return None


def verify_all_seven_fish(verbose=True, output_path="results_cascade_n7.txt"):
    """
    For each of the seven n=7 fish multisets, compute its fingerprint
    equation directly, then verify each cousin dies by layer-0 (or, where
    needed, by a one-step Laurent kill at a neighboring zone).

    Writes the cascade trace to `output_path`.
    """
    n = N_VALUE
    log_lines = []

    def log(s):
        log_lines.append(s)
        if verbose:
            print(s)

    log("=" * 78)
    log("Cascade verification of the seven n=7 'fish' multisets")
    log("=" * 78)
    log("")
    log("For each fish multiset M (a non-triangulation chord multiset of")
    log("the heptagon not killed by the layer-0 mechanism at any zone),")
    log("we exhibit a zone Z and Laurent order k such that the fingerprint")
    log("equation at (Z, k) -- after substituting in lower-tier kills --")
    log("collapses to a_M = 0.  All Laurent expansions are computed")
    log("symbolically with sympy.")
    log("")

    all_dead = True

    for i, fish in enumerate(SEVEN_FISH):
        log("-" * 78)
        log(f"Fish #{i+1}: M = {format_multiset(fish)}")
        log("-" * 78)
        rec = verify_each_fish_from_scratch(fish, n, verbose=False)
        if rec is None:
            log(f"  !!! No kill recipe found for fish #{i+1}.")
            all_dead = False
            continue
        zone_r = rec["zone_r"]
        sp_z, _, _ = zone_structure(zone_r, n)
        log(f"  Kill zone Z = Z_{{{sp_z[0]},{sp_z[1]}}};  Laurent order k = {rec['order']};")
        # Pretty fingerprint
        fp = rec["fingerprint"]
        log(f"  fingerprint = {sp.simplify(fp)}.")
        log(f"  (cyclic shift v -> v + {rec['shift']} (mod {n}) maps fish #1 to this fish.)")
        log("")
        log(f"  At this zone+order+fingerprint, the equation reads:")
        log("")
        for ms, scalar in rec["equation"]:
            marker = "  <-- fish" if tuple(sorted(ms)) == tuple(sorted(fish)) else ""
            sign = "+" if scalar >= 0 else ""
            log(f"     {sign}{scalar}  *  a_{format_multiset(ms)}{marker}")
        log("")
        # Verify every cousin dies.
        cousins = [(m, s) for m, s in rec["equation"]
                   if tuple(sorted(m)) != tuple(sorted(fish))]
        all_cousins_dead = True
        log(f"  Cousins (must die before this equation can isolate the fish):")
        for ms, scalar in cousins:
            z0 = layer0_kill_zone(ms, n)
            if z0 is not None:
                sp_c, _, _ = zone_structure(z0, n)
                log(f"     a_{format_multiset(ms)}  =  0  by LAYER-0 at "
                    f"Z_{{{sp_c[0]},{sp_c[1]}}}.")
            else:
                # Try Laurent kill at any zone.
                killed_by_laurent = False
                for try_r in range(1, n + 1):
                    rl = verify_at_zone_laurent_kill(ms, n, try_r,
                                                    target_order=2, verbose=False)
                    if rl is not None:
                        sp_c, _, _ = zone_structure(try_r, n)
                        comp = rl["companion"]
                        log(f"     a_{format_multiset(ms)}  =  0  by Z_{{{sp_c[0]},{sp_c[1]}}} "
                            f"Laurent kill (order X_{sp_c[0]}{sp_c[1]}^-2, "
                            f"fingerprint X_{comp[0]}{comp[1]}/(...)).")
                        killed_by_laurent = True
                        break
                if not killed_by_laurent:
                    log(f"     a_{format_multiset(ms)}  :  !!! UNVERIFIED !!!")
                    all_cousins_dead = False
        log("")
        # Conclusion.
        fish_scalar = next(s for m, s in rec["equation"]
                           if tuple(sorted(m)) == tuple(sorted(fish)))
        log(f"  Fish term in equation: ({fish_scalar}) * a_M.")
        if all_cousins_dead:
            log(f"  All cousins killed -> equation collapses to "
                f"({fish_scalar}) * a_M = 0  =>  a_M = 0.")
            log(f"  QED for fish #{i+1}.")
        else:
            log(f"  Not all cousins killed -- this fish's verification is INCOMPLETE.")
            all_dead = False
        log("")

    log("=" * 78)
    if all_dead:
        log("Summary: all seven n=7 fish coefficients are killed by the cascade.")
        log("Layer-0 alone does NOT suffice for any of them; the Laurent tail")
        log("at a neighboring zone provides the additional fingerprint that")
        log("isolates each fish once its cousins die at lower tier.")
    else:
        log("Summary: at least one fish was not fully verified.")
    log("=" * 78)

    with open(output_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\nResults written to {output_path}")
    return log_lines


def verify_z46_laurent_kill_rotated(ms_rot, n, shift, verbose=False):
    """Verify the cyclic rotation of the Z_{4,6} Laurent kill."""
    zone_r_rot = ((4 - 1 + shift) % n) + 1
    return verify_at_zone_laurent_kill(ms_rot, n, zone_r_rot, verbose=verbose)


def verify_at_zone_laurent_kill(ms, n, zone_r, target_order=2, verbose=False):
    """Generalized: try Laurent-tail kill at zone Z_{zone_r,...} at given order."""
    special, bare, pairs = zone_structure(zone_r, n)
    substitutes_in_ms = []
    for c in ms:
        for companion, substitute in pairs:
            if c == substitute:
                substitutes_in_ms.append((companion, substitute))
                break
    if not substitutes_in_ms:
        return None
    for companion, substitute in substitutes_in_ms:
        ms_minus = list(ms)
        ms_minus.remove(substitute)
        rest_factors = sp.Integer(1)
        for c in ms_minus:
            rest_factors *= chord_sym(c)
        fp = chord_sym(companion) / rest_factors
        pool = candidate_multisets_for_fingerprint(fp, n, n - 3)
        eq = fingerprint_equation(ms, zone_r, n, target_order, fp, pool)
        contributors = [(m, s) for m, s in eq if s != 0]
        if len(contributors) == 1 and contributors[0][0] == ms:
            return {"zone_r": zone_r, "order": target_order,
                    "companion": companion, "equation": eq}
    return None


def format_multiset(ms):
    """Pretty-print a multiset as {(i,j),(k,l),...}."""
    return "{" + ",".join(f"({i},{j})" for i, j in ms) + "}"


if __name__ == "__main__":
    verify_all_seven_fish(verbose=True, output_path="results_cascade_n7.txt")
