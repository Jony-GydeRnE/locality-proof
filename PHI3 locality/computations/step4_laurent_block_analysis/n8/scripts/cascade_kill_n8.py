"""
================================================================================
cascade_kill_n8.py  --  Search for Laurent-cascade kills of every n=8 step-1
                        survivor (the 100 "fish" multisets at n=8).
================================================================================

PURPOSE (BIG PICTURE)
---------------------
At n=8, the leading-order ("layer-0" / "Step-1") hidden-zero kill mechanism
fails to kill exactly 100 non-triangulation chord multisets of size n-3 = 5.
This script asks: does the same kind of LAURENT CASCADE that worked at n=7
(file ../cascade_n7/cascade_kill_n7.py) also kill every n=8 survivor?

At n=7 we had a single hand-picked recipe for "fish #1" that, by cyclic
rotation, gave the recipe for the other six fish. At n=8 we cannot rely on
a hand-picked recipe -- there are 100 survivors falling into multiple cyclic
orbits, and a priori we do not know which (zone, Laurent order, fingerprint)
choice will isolate each survivor.  So instead of a hand-picked recipe, this
script SEARCHES for one per survivor, trying:

  for each survivor M:
    for each zone Z = Z_{r,r+2}:
      if M has at least one (companion, substitute) pair structure on Z:
        for each substitute c_sub in M (with companion c_comp):
          for each Laurent order k in {leading_order + 1, leading_order + 2}:
            build fingerprint  fp = X_{c_comp} / (product of the rest of M's chords)
            compute the fingerprint equation at (Z, k, fp);
            if M appears with nonzero coefficient AND every cousin in the
              equation is step-1-killable (or itself a survivor with a
              previously-found cascade recipe), record this as a kill recipe.

If a depth-1 cascade is found for every survivor, the conclusion is:

   *Every n=8 step-1 survivor coefficient is zero, by a depth-1 Laurent
    cascade with prerequisites all of step-1 type.*

If some survivors do NOT admit a depth-1 cascade, the script reports them
explicitly: they are the candidates for needing a deeper cascade (or for
revealing that the depth-1 phenomenon stops at n=7).

CODE STYLE NOTE
---------------
Per the repo contribution guidelines, every function below has a docstring
that explains both the CODE LOGIC (what algorithm is implemented) and the
PHYSICS / MATHEMATICS (what the algorithm computes in the language of the
proof).  Math/physics readers should be able to follow the logic without
ever needing to read the code line-by-line; programmers should be able to
follow the implementation without ever needing to read the proof.

SHARED UTILITIES
----------------
The combinatorial primitives (chord normalization, zone structure, layer-0
killability) are repeated verbatim from cascade_kill_n7.py so this script
is self-contained.  If you change conventions in one place, change them in
both.
"""

from itertools import combinations_with_replacement
import sympy as sp
import time
import sys


# ============================================================================
# 1.  COMBINATORICS:  chords, zones, layer-0 killability
#     (verbatim from cascade_kill_n7.py)
# ============================================================================

def normalize(i, j, n):
    """
    Canonical name for the chord {i, j} of the n-gon.

    LOGIC:  reduce both endpoints into {1, ..., n} mod n, then sort so the
            smaller one comes first.  Returns a tuple (a, b) with a < b.

    PHYSICS:  X_{ij} = X_{ji} kinematically, so we always store as the
              ordered pair (a, b) with a < b for hashing/equality.
    """
    i = ((i - 1) % n) + 1
    j = ((j - 1) % n) + 1
    if i > j:
        i, j = j, i
    return (i, j)


def all_chords(n):
    """
    All planar chords of the n-gon (non-edge pairs).

    LOGIC:  iterate over pairs (i, j) with 1 <= i < j <= n, drop the cases
            where the pair is a polygon edge (j - i == 1) or the wrap-around
            edge (i, j) == (1, n).

    PHYSICS:  these are the d = n(n-3)/2 planar Mandelstam invariants
              X_{ij} = (p_i + ... + p_{j-1})^2 of the colour-ordered amplitude.
    """
    return [(i, j) for i in range(1, n + 1)
                   for j in range(i + 1, n + 1)
                   if j - i >= 2 and not (i == 1 and j == n)]


def zone_structure(r, n):
    """
    For zone Z_{r, r+2}, return (special, bare, pairs).

    LOGIC:  identify the special chord (r, r+2), the bare chord (r+1, r-1)
            (with r-1 taken cyclically), and the n-4 (companion, substitute)
            pairs (X_{r,k}, X_{r+1,k}) for k = r+3, ..., r+n-2 (cyclic).

    PHYSICS:  on the codimension-2 hidden-zero locus where all c_{r,*}
              vanish, the master substitution rewrites every chord
              incident to vertex r+1 as a Q-linear function of free
              chords + the special chord:
                  X_{r+1, k}    = X_{r, k} - X_{r, r+2}     (substitute)
                  X_{r+1, r-1}  = -X_{r, r+2}               (bare)
              The "companion" of the substitute X_{r+1,k} is the row-r
              chord X_{r,k}; both have endpoint r in common.
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
    Layer-0 / Step-1 killability of `ms_list` at the given zone.

    LOGIC:  return True iff
              (K1) the multiset contains neither the special nor the bare;
              (K2) for every (companion, substitute) pair, the multiset
                   does NOT contain BOTH chords of the pair.

    PHYSICS:  if both (K1) and (K2) hold, then sending the special chord
              to infinity at this zone forces the coefficient a_M to be
              zero by the leading-Laurent argument (linear independence
              of monomials in the free variables).
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
    """
    Smallest zone index r such that `ms_list` is layer-0 killable at
    Z_{r, r+2}, or None if no zone works.

    LOGIC:  scan r = 1, 2, ..., n; return the first hit.

    PHYSICS:  if any zone kills the multiset, we say it is "step-1
              killable"; the survivors are exactly the multisets that
              fail (K1) or (K2) at every cyclic zone.
    """
    for r in range(1, n + 1):
        if killable_at_zone_layer0(ms_list, zone_structure(r, n)):
            return r
    return None


def is_triangulation(ms_list, n):
    """
    Test whether the multiset is a triangulation of the n-gon.

    LOGIC:  triangulations of an n-gon have exactly n-3 distinct,
            pairwise non-crossing chords.  Two chords (a, b) and (c, d)
            cross iff (in cyclic order around the polygon) one of c, d
            is strictly between a and b and the other is not.  We check
            distinctness (multiset = set of size n-3) and pairwise
            non-crossing.

    PHYSICS:  triangulation coefficients in B are exactly the Feynman
              diagram coefficients of A_n^tree.  We exclude them from the
              survivor pool because Theorem 1 of triangulation_layer0.tex
              shows they are NEVER step-1 killable, so they never appear
              in a step-1 survivor list anyway -- but the explicit check
              is a sanity guard.
    """
    if len(set(ms_list)) != n - 3:
        return False  # repeated chord => double pole, not a triangulation
    chords = list(set(ms_list))
    for i in range(len(chords)):
        for j in range(i + 1, len(chords)):
            a, b = chords[i]
            c, d = chords[j]
            # If the chords share an endpoint, they cannot cross.
            # (Crossing requires four distinct vertices.)
            if a == c or a == d or b == c or b == d:
                continue
            # Crossing test: exactly one of c, d strictly between a and b.
            # With distinct endpoints and a<b, c<d, this is the standard
            # convex-polygon chord-crossing condition.
            between_c = (a < c < b)
            between_d = (a < d < b)
            if between_c != between_d:
                return False
    return True


# ============================================================================
# 2.  SYMBOLIC LAURENT EXPANSION OF A MONOMIAL ON A ZONE
#     (verbatim from cascade_kill_n7.py)
# ============================================================================

def chord_symbol(c):
    """
    Sympy positive symbol X_{i,j} for chord c = (i, j).

    LOGIC:  return sp.Symbol(f"X_{i}_{j}", positive=True).

    PHYSICS:  X_{ij} > 0 in the kinematic region where these polynomial
              identities are most easily checked; the positive=True hint
              also helps sympy simplify Laurent expansions correctly.
    """
    i, j = c
    return sp.Symbol(f"X_{i}_{j}", positive=True)


def laurent_coefficient(monomial_chords, zone_r, n, max_order):
    """
    Laurent expansion of  1 / prod_{c in monomial_chords} X_c  on zone
    Z_{zone_r, zone_r+2}, up to order X_S^{-max_order} where X_S is the
    special chord of the zone.

    LOGIC:
      1. Substitute every BARE chord by  -X_S.
      2. Substitute every SUBSTITUTE chord X_{r+1,k} by  (X_{r,k} - X_S).
      3. Free chords (companions and born-free) stay as-is.
      4. Replace X_S = 1/u and use sp.series(..., u, 0, max_order+1) to
         expand the resulting expression around u = 0.
      5. Read off coefficients of u^k = X_S^{-k} for k = 0, ..., max_order.

    PHYSICS:  this is the Laurent expansion of the rational function on
              the codim-2 hidden-zero locus, expanded in the small
              parameter 1/X_S.  Algebraic independence of the free
              variables means each coefficient at each Laurent order is
              an independent constraint on the ansatz coefficients a_M.

    Returns a dict  { order_k : sympy expression in free variables }.
    """
    special, bare, pairs = zone_structure(zone_r, n)
    X_S = chord_symbol(special)
    u   = sp.Symbol("__u__", positive=True)

    # Step 1-3: build the substituted monomial.
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
                expr = expr / chord_symbol(c)

    # Step 4-5: expand in 1/X_S via the substitution X_S = 1/u.
    expr_u = expr.subs(X_S, 1 / u)
    s = sp.series(expr_u, u, 0, max_order + 1).removeO()
    s = sp.expand(s)
    out = {}
    for k in range(max_order + 1):
        coeff = sp.simplify(s.coeff(u, k))
        out[k] = coeff
    return out


def coefficient_of_monomial(expr, target_monomial, free_vars):
    """
    Extract the scalar coefficient of `target_monomial` from `expr`,
    where both are rational expressions in `free_vars`.

    LOGIC:  divide expr by target_monomial, expand into a sum, collect
            terms that are constant (no free-variable dependence).

    PHYSICS:  picks the integer/rational scalar in front of one specific
              free-variable monomial (the "fingerprint" of an ansatz
              monomial at a chosen Laurent order).
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


def fingerprint_equation(target_multiset, zone_r, n, target_order,
                         target_fingerprint, candidates_pool):
    """
    Compute the equation  sum_{M' in pool} scalar(M') * a_{M'} = 0  where
    scalar(M') is the coefficient of `target_fingerprint` in M''s Laurent
    expansion at order `target_order` on Z_{zone_r, zone_r+2}.

    LOGIC:  iterate over `candidates_pool`, expand each via
            laurent_coefficient, extract the fingerprint coefficient via
            coefficient_of_monomial, accumulate (multiset, scalar) pairs
            with nonzero scalar.

    PHYSICS:  the equation is the Laurent-tail constraint that
              B|_{Z_{r,r+2}} = 0 at order X_S^{-target_order} on the
              specific free-variable monomial `target_fingerprint`.  Any
              ansatz monomial that contributes to that order on that
              fingerprint enters the equation with the integer scalar
              produced by the geometric series 1/(X_companion - X_S) =
              -1/X_S - X_companion/X_S^2 - ...

    Returns a list of (multiset, scalar) pairs (each scalar nonzero).
    """
    equation_terms = []
    for ms in candidates_pool:
        coeff_dict = laurent_coefficient(ms, zone_r, n, target_order)
        order_k_expr = coeff_dict[target_order]
        if order_k_expr == 0:
            continue
        # Compute free vars on this zone (companions + born-free).
        special, bare, pairs = zone_structure(zone_r, n)
        free_vars = []
        for c in all_chords(n):
            if c == special:
                continue
            sym = chord_symbol(c)
            is_substitute = any(c == sub for _, sub in pairs)
            is_bare       = (c == bare)
            if is_substitute or is_bare:
                continue
            free_vars.append(sym)
        scalar = coefficient_of_monomial(order_k_expr, target_fingerprint,
                                         free_vars)
        if scalar is None or scalar == 0:
            continue
        equation_terms.append((ms, scalar))
    return equation_terms


def candidate_multisets_for_fingerprint(target_fp, n, multiset_size):
    """
    Enumerate the multisets that COULD contribute the rational monomial
    `target_fp` at SOME Laurent order on SOME zone.

    LOGIC:  read off the denominator chords of `target_fp` (with
            multiplicities) -- these MUST appear in any contributing
            multiset.  Add any combinatorial choice of extra chords to
            reach `multiset_size`.

    PHYSICS:  a fast pre-filter that drops the candidate-pool size from
              C(d + N - 1, N) (full ansatz) to roughly d^k (where k is
              the number of "extra slots" beyond the forced denominator),
              making the search tractable at n = 8.  False positives (a
              candidate that has the right denominator chords but does
              not actually produce `target_fp` at any zone) are pruned
              later by coefficient_of_monomial returning 0.
    """
    fp_clean = sp.simplify(target_fp)
    num, den = sp.fraction(fp_clean)
    den_syms_count = {}
    if den != 1:
        for term in sp.Mul.make_args(den):
            base, exp = term.as_base_exp()
            if base.is_Symbol:
                den_syms_count[base] = den_syms_count.get(base, 0) + int(exp)

    def sym_to_chord(s):
        name = s.name
        parts = name.split("_")
        return (int(parts[1]), int(parts[2]))

    den_chord_count = {sym_to_chord(s): k for s, k in den_syms_count.items()}
    forced_chords = []
    for c, k in den_chord_count.items():
        forced_chords.extend([c] * k)
    extras_needed = multiset_size - len(forced_chords)
    if extras_needed < 0:
        return []

    chords = all_chords(n)
    seen = set()
    deduped = []
    for extras in combinations_with_replacement(chords, extras_needed):
        candidate = tuple(sorted(forced_chords + list(extras)))
        if candidate not in seen:
            seen.add(candidate)
            deduped.append(candidate)
    return deduped


# ============================================================================
# 3.  ENUMERATE THE n=8 STEP-1 SURVIVORS
# ============================================================================

def all_size_N_multisets(n):
    """
    Every chord multiset of size n-3.

    LOGIC:  combinations_with_replacement of all_chords(n) taken n-3 at
            a time.

    PHYSICS:  these index every monomial of the most general planar
              rational ansatz B with at most simple poles at mass
              dimension -2(n-3).  Sizes:
                 n = 5:    15
                 n = 6:   165
                 n = 7:  2380
                 n = 8: 42504
                 n = 9: 906192
    """
    chords = all_chords(n)
    return list(combinations_with_replacement(chords, n - 3))


def step1_survivors(n):
    """
    Compute the list of step-1 SURVIVORS at n: multisets that fail the
    layer-0 kill at every zone, AND that are not triangulations.

    LOGIC:  full enumeration; for each multiset, call layer0_kill_zone;
            keep those where it returns None and which are NOT
            triangulations.

    PHYSICS:  these are the "fish" -- the input set this script must
              kill via the Laurent cascade.  At n = 8, there should be
              exactly 100 of them (per the kill_enumeration headline
              numbers in computations/step1_layer0_kill/).
    """
    out = []
    for ms in all_size_N_multisets(n):
        if layer0_kill_zone(list(ms), n) is None:
            if not is_triangulation(list(ms), n):
                out.append(ms)
    return out


# ============================================================================
# 4.  RECIPE SEARCH: per-fish depth-1 cascade
# ============================================================================
#
# Strategy (depth-1 cascade): for each survivor M, we look for a triple
# (zone Z, Laurent order k, fingerprint fp) such that:
#
#   (a) M appears in the fingerprint equation at (Z, k, fp) with nonzero
#       integer scalar;
#   (b) every OTHER multiset M' in the equation is layer-0 killable
#       (at some zone -- not necessarily the same as Z).
#
# Then the equation, after substituting the cousin kills, collapses to
# scalar(M) * a_M = 0, hence a_M = 0.
#
# We generate fingerprint candidates by mimicking the n=7 worked example:
# for each (companion, substitute) pair (X_comp, X_sub) with X_sub in M,
# the fingerprint  fp = X_comp / (M - X_sub)  is the natural one produced
# by the Laurent tail of 1/(X_comp - X_S).  We try the leading-tail order
# (#bares + #substitutes_in_M + 1) and one order deeper.

def fish_substitutes_at_zone(fish, zone_r, n):
    """
    Identify which chords of `fish` are substitutes at zone Z_{zone_r,...}.

    LOGIC:  for each chord c in fish, scan the (companion, substitute)
            pairs of the zone; if c equals a substitute, record the pair.

    PHYSICS:  substitutes are exactly the chords with vertex r+1 as
              endpoint, except the bare.  Their Laurent expansion via
              1/(X_comp - X_S) = -sum_{k>=1} X_comp^{k-1} / X_S^k is
              what produces all the nontrivial "fingerprints" beyond
              leading order.
    """
    _, _, pairs = zone_structure(zone_r, n)
    out = []
    for c in fish:
        for companion, substitute in pairs:
            if c == substitute:
                out.append((companion, substitute))
                break
    return out


def fish_count_bare_and_subs(fish, zone_r, n):
    """
    Count occurrences (with multiplicity) of the bare chord and of any
    substitute chord in `fish`, on zone Z_{zone_r, zone_r+2}.

    LOGIC:  sum 1 for each occurrence in the multiset.

    PHYSICS:  the leading Laurent order in 1/X_S of the fish's monomial
              on this zone is exactly  (# bares) + (# substitutes), since
              each bare and each substitute contributes one factor of
              1/X_S to the leading expansion.  So the first NONTRIVIAL
              Laurent order (the one with companion variables in it) is
              the next one up.
    """
    special, bare, pairs = zone_structure(zone_r, n)
    n_bare = sum(1 for c in fish if c == bare)
    n_sub  = sum(1 for c in fish if any(c == sub for _, sub in pairs))
    return n_bare, n_sub


def search_depth1_recipe(fish, n, max_extra_orders=2, verbose=False):
    """
    Search for a depth-1 cascade kill recipe for `fish`.

    LOGIC:
      For zone_r in 1..n:
        compute leading_order = (#bares + #subs) of fish on this zone.
        For each substitute (c_comp, c_sub) in fish on this zone:
          For extra in 1..max_extra_orders:
            target_order = leading_order + extra
            build fingerprint  fp = X_{c_comp} / (X factors of fish minus c_sub)
            generate candidates_pool = candidate_multisets_for_fingerprint(fp)
            equation = fingerprint_equation(fish, zone_r, n, target_order, fp,
                                             candidates_pool)
            if fish in equation with nonzero scalar AND every other
              multiset in equation is layer-0 killable somewhere:
                return this recipe.
      return None.

    PHYSICS:  exactly the n=7 cascade structure ("one Laurent step deeper
              at one zone, with prerequisites all of layer-0 type at one
              neighbouring zone"), generalized to n=8 with no a-priori
              hand-picked recipe.

    Returns either None or a dict with keys:
      'zone_r', 'order', 'fingerprint', 'companion', 'equation',
      'cousin_kills' = list of (multiset, scalar, kill_zone).
    """
    for zone_r in range(1, n + 1):
        n_bare, n_sub = fish_count_bare_and_subs(fish, zone_r, n)
        if n_sub == 0:
            # No substitutes at this zone -> Laurent tail produces no companion
            # variables; skip.
            continue
        leading_order = n_bare + n_sub
        substitutes_in_fish = fish_substitutes_at_zone(fish, zone_r, n)
        # Deduplicate substitutes by their (companion, substitute) identity:
        # if the same substitute appears with multiplicity, the fingerprint
        # X_comp / (rest minus one occurrence) is the same.
        seen_pairs = set()
        unique_pairs = []
        for pair in substitutes_in_fish:
            if pair not in seen_pairs:
                seen_pairs.add(pair)
                unique_pairs.append(pair)
        # On this zone, identify the special, bare, and substitute chord set
        # so we can exclude them from the fingerprint denominator (those get
        # replaced by powers of X_S, not free-variable factors).
        special_z, bare_z, pairs_z = zone_structure(zone_r, n)
        substitute_set = {sub for _, sub in pairs_z}
        for companion, substitute in unique_pairs:
            # Build the fingerprint by:
            #   numerator   = X_{companion of the chosen substitute}
            #   denominator = product of X_c for c a FREE chord of the fish,
            #                 minus the chosen substitute (and minus all
            #                 bares, the special, and other substitutes,
            #                 since those don't contribute as free vars at
            #                 order leading_order + 1 -- they contribute
            #                 only powers of X_S).
            #
            # We remove ONE occurrence of the chosen substitute from the
            # multiset, then drop any chord that is the special, the bare,
            # or any other substitute on this zone.
            ms_minus = list(fish)
            ms_minus.remove(substitute)
            rest_free = [c for c in ms_minus
                         if c != special_z
                         and c != bare_z
                         and c not in substitute_set]
            rest_factors = sp.Integer(1)
            for c in rest_free:
                rest_factors *= chord_symbol(c)
            fp = chord_symbol(companion) / rest_factors
            for extra in range(1, max_extra_orders + 1):
                target_order = leading_order + extra
                pool = candidate_multisets_for_fingerprint(fp, n, n - 3)
                if not pool:
                    continue
                eq = fingerprint_equation(fish, zone_r, n, target_order, fp, pool)
                if not eq:
                    continue
                # Find fish in equation.
                fish_scalar = None
                fish_sorted = tuple(sorted(fish))
                for m, s in eq:
                    if tuple(sorted(m)) == fish_sorted:
                        fish_scalar = s
                        break
                if fish_scalar is None or fish_scalar == 0:
                    continue
                # Check every cousin is layer-0 killable.
                cousins = [(m, s) for m, s in eq
                           if tuple(sorted(m)) != fish_sorted]
                cousin_kills = []
                all_dead = True
                for m, s in cousins:
                    z0 = layer0_kill_zone(list(m), n)
                    if z0 is None:
                        all_dead = False
                        cousin_kills.append((m, s, None))
                    else:
                        cousin_kills.append((m, s, z0))
                if all_dead:
                    if verbose:
                        print(f"   [hit] zone Z_{{{zone_r},{(zone_r+2-1)%n+1}}}, "
                              f"order {target_order}, "
                              f"companion {companion}: equation has "
                              f"{len(eq)} terms, all cousins die by step-1.")
                    return {
                        'zone_r': zone_r,
                        'order': target_order,
                        'fingerprint': fp,
                        'companion': companion,
                        'fish_scalar': fish_scalar,
                        'equation': eq,
                        'cousin_kills': cousin_kills,
                    }
    return None


# ============================================================================
# 5.  MAIN DRIVER
# ============================================================================

def format_multiset(ms):
    """Pretty-print a multiset as {(i,j),(k,l),...}."""
    return "{" + ",".join(f"({i},{j})" for i, j in ms) + "}"


def main(output_path="results_cascade_n8.txt", limit=None):
    """
    Top-level driver.  Compute n=8 step-1 survivors, search for a depth-1
    cascade recipe per survivor, write a human-readable trace.

    LOGIC:
      1. Enumerate the 100 step-1 survivors of size n-3 = 5.
      2. For each survivor, call search_depth1_recipe.
      3. Record SUCCESS (recipe found, all cousins step-1 killable) or
         FAILURE (no depth-1 recipe found).
      4. Print summary: # survivors killed by depth-1, # remaining.

    PHYSICS:  if the success rate is 100%, this confirms the conjecture
              (from cascade_n7.tex's outlook) that every step-1 survivor
              at n >= 7 dies via a depth-1 Laurent cascade with all
              prerequisites of layer-0 type.  If success rate < 100%,
              the failures are the candidates needing a deeper cascade,
              and the script's job is to identify them precisely.
    """
    n = 8
    log_lines = []

    def log(s):
        log_lines.append(s)
        print(s)
        sys.stdout.flush()

    log("=" * 78)
    log(f"Cascade verification of the n=8 step-1 survivors")
    log("=" * 78)
    log("")
    log(f"Strategy: for each survivor M, search for a depth-1 cascade recipe.")
    log(f"A recipe is (zone Z, Laurent order k, fingerprint fp) such that:")
    log(f"  (a) M appears in the fingerprint equation at (Z, k, fp) with")
    log(f"      nonzero integer scalar;")
    log(f"  (b) every other multiset in the equation is step-1 killable")
    log(f"      at some zone (not necessarily Z).")
    log(f"If such a recipe is found, a_M = 0 follows after substitution.")
    log("")

    log("Enumerating size-5 multisets and computing step-1 survivors...")
    t0 = time.time()
    survivors = step1_survivors(n)
    log(f"  Found {len(survivors)} step-1 survivors in "
        f"{time.time()-t0:.1f} s (expected: 100).")
    log("")

    if limit is not None:
        survivors = survivors[:limit]
        log(f"NOTE: limit set; only checking the first {len(survivors)} "
            f"survivors.")
        log("")

    successes = []
    failures = []
    t_start = time.time()
    for i, fish in enumerate(survivors):
        log("-" * 78)
        log(f"Survivor #{i+1}/{len(survivors)}: M = {format_multiset(fish)}")
        log("-" * 78)
        t_fish = time.time()
        rec = search_depth1_recipe(fish, n, max_extra_orders=2, verbose=False)
        elapsed = time.time() - t_fish
        if rec is None:
            log(f"  !!! No depth-1 recipe found  ({elapsed:.1f} s).")
            failures.append(fish)
            log("")
            continue
        zone_r = rec['zone_r']
        sp_z, _, _ = zone_structure(zone_r, n)
        log(f"  Kill zone Z = Z_{{{sp_z[0]},{sp_z[1]}}};  Laurent order k = "
            f"{rec['order']};")
        log(f"  fingerprint = {sp.simplify(rec['fingerprint'])}  "
            f"(companion {rec['companion']}).")
        log(f"  ({elapsed:.1f} s)")
        log("")
        log("  Equation:")
        for ms, s in rec['equation']:
            marker = "  <-- fish" if tuple(sorted(ms)) == tuple(sorted(fish)) else ""
            sign = "+" if s >= 0 else ""
            log(f"     {sign}{s}  *  a_{format_multiset(ms)}{marker}")
        log("")
        log("  Cousin step-1 kills (all required for collapse to a_M = 0):")
        for ms, s, z0 in rec['cousin_kills']:
            sp_c, _, _ = zone_structure(z0, n)
            log(f"     a_{format_multiset(ms)}  =  0  by step-1 at "
                f"Z_{{{sp_c[0]},{sp_c[1]}}}.")
        log("")
        log(f"  Fish term in equation: ({rec['fish_scalar']}) * a_M.")
        log(f"  All cousins killed -> equation collapses to "
            f"({rec['fish_scalar']}) * a_M = 0  =>  a_M = 0.")
        log(f"  QED for survivor #{i+1}.")
        log("")
        successes.append((fish, rec))

    total_elapsed = time.time() - t_start
    log("=" * 78)
    log(f"Summary: {len(successes)} / {len(survivors)} survivors killed by "
        f"depth-1 cascade in {total_elapsed:.1f} s.")
    if failures:
        log(f"Failures ({len(failures)}):")
        for f in failures:
            log(f"  {format_multiset(f)}")
    else:
        log("All n=8 step-1 survivors die via depth-1 Laurent cascade.")
        log("Conclusion: every n=8 non-triangulation coefficient a_M is")
        log("forced to zero by step-1 (where applicable) or by step-1 plus")
        log("a single Laurent-tail step at one neighbouring zone.")
    log("=" * 78)

    with open(output_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    print(f"\nResults written to {output_path}")
    return successes, failures


if __name__ == "__main__":
    # Optional CLI: pass an integer to limit the number of survivors checked
    # (useful for quick smoke tests).
    limit = None
    if len(sys.argv) > 1:
        limit = int(sys.argv[1])
    main(output_path="results_cascade_n8.txt", limit=limit)
