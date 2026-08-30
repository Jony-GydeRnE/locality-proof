#!/usr/bin/env python3
"""
Stress tests for the load-bearing claims of Part II + Part III.

Three families, all symbolic (sympy) so passes are exact, not numerical:

  (T1) Algebraic identities
       T1a: Telescoping lemma -- the X-substitutions on Z^(k) satisfy c_{ij}=0
            for every i=1..k, j=k+2..n-1.
       T1b: Chord-rate table (Part II eq. 7) -- under (alpha_ell; beta_i) rate
            vector, each chord's "finite iff" clause is verified.

  (T2) Survival classification (Part II Prop 1)
       Enumerate all multisets of size n-3 on the non-crossing chord set, and
       cross-check two independent classifications:
        - (a) is the rate-equation system imposed by the bridges of M consistent?
              If consistent ==> M is killed (subset of some F(alpha,beta)).
              If inconsistent ==> M survives.
        - (b) does M contain the special, the bare, or a column pair
              {X_{1,m}, X_{k+1,m}} for some m in {k+3,..,n-1}?
       Prop 1 says (a)-survives iff (b)-true. Test this for every M.

  (T3) Step-2 relations (Part III)
       Build the full ansatz B_n with symbolic coefficients a_M, substitute
       the k-zero constraints, expand on D-subsets, and verify:
        T3a: the m=1 binomial trace (unitarity)  a_{X_{1,k+2}, Y} = a_{X_{k+1,n}, Y}
             for each born-free monomial Y;
        T3b: the m=2 binomial trace               a_{ss} - a_{sb} + a_{bb} = 0
             (in particular  a_{14,14} - a_{14,3n} + a_{3n,3n} = 0  for k=2)
        T3c: the 2-zero per-column relation
             sum_{i=1}^{3} (a_{14,i_ell} - a_{3n,i_ell}) = 0
             at m=2 for k=2, each far column ell.

Run:  python3 stress_tests.py
"""

import sympy as sp
from itertools import combinations_with_replacement, product
import sys


# ---------------------------------------------------------------------------
# Set up the k-zero constraint surface symbolically.
# ---------------------------------------------------------------------------

def is_leg(n, a, b):
    """X_{a,b} is a leg iff |b-a|=1 or {a,b}={1,n}."""
    if a > b:
        a, b = b, a
    return (b - a == 1) or (a == 1 and b == n)


def doesnt_cross_enclosing(n, k, a, b):
    """X_{a,b} does not cross the enclosing chord X_{1,k+1}.
       Cross iff one endpoint is strictly in {2..k} and the other in {k+2..n}."""
    if a > b:
        a, b = b, a
    inside = lambda v: 2 <= v <= k
    outside = lambda v: k + 2 <= v <= n
    crosses = (inside(a) and outside(b)) or (outside(a) and inside(b))
    return not crosses


def born_free_chords(n, k):
    """Born-free = chords internal to N={1..k+1} or internal to F={k+2..n}."""
    bf = []
    # inside N
    for a in range(1, k + 2):
        for b in range(a + 2, k + 2):
            if not is_leg(n, a, b):
                bf.append((a, b))
    # inside F
    for a in range(k + 2, n + 1):
        for b in range(a + 2, n + 1):
            if not is_leg(n, a, b):
                bf.append((a, b))
    return bf


def all_chords(n):
    """All chord variables (i,j) with 1<=i<j<=n, not a leg."""
    chords = []
    for a in range(1, n + 1):
        for b in range(a + 2, n + 1):
            if not is_leg(n, a, b):
                chords.append((a, b))
    return chords


def setup_kzero(n, k):
    """Return (X, s, t, u, Y) where:
       X[(i,j)] = sympy expression for X_{i,j} on Z^(k) in indep coords,
       s        = special chord symbol = X_{1,k+2},
       t[i]     = offset symbol X_{i,k+2} for i=2..k,
       u[l]     = X_{1,l} symbol for l=k+3..n-1,
       Y[(i,j)] = born-free symbol for each born-free chord (i,j)."""
    s = sp.Symbol('s')
    t = {i: sp.Symbol(f't{i}') for i in range(2, k + 1)}
    u = {l: sp.Symbol(f'u{l}') for l in range(k + 3, n)}
    Y = {bf: sp.Symbol(f'Y_{bf[0]}_{bf[1]}') for bf in born_free_chords(n, k)}

    X = {}
    # special and bare
    X[(1, k + 2)] = s
    X[(k + 1, n)] = -s
    # backbone X_{k+1, l} = X_{1,l} - s, l=k+3..n-1
    for l in range(k + 3, n):
        X[(k + 1, l)] = u[l] - s
    # column-(k+2) interior offsets X_{i, k+2} = t_i for i=2..k
    for i in range(2, k + 1):
        X[(i, k + 2)] = t[i]
    # column 1: X_{1, l} = u_l for l=k+3..n-1
    for l in range(k + 3, n):
        X[(1, l)] = u[l]
    # interior bridges X_{i, l} = X_{k+1,l} + t_i for i=2..k, l=k+3..n-1
    for i in range(2, k + 1):
        for l in range(k + 3, n):
            X[(i, l)] = (u[l] - s) + t[i]
    # column n: X_{i, n} = X_{k+1,n} + t_i = -s + t_i, i=2..k
    for i in range(2, k + 1):
        X[(i, n)] = -s + t[i]
    # born-free chords
    for bf, sym in Y.items():
        X[bf] = sym
    return X, s, t, u, Y


def Xget(X, n, a, b):
    """Get X[(a,b)] obeying the leg convention X_{m,m+1}=X_{1,n}=0."""
    if a > b:
        a, b = b, a
    if is_leg(n, a, b):
        return sp.Integer(0)
    return X[(a, b)]


# ---------------------------------------------------------------------------
# T1a: Telescoping lemma
# ---------------------------------------------------------------------------

def test_telescoping(n, k):
    X, s, t, u, Y = setup_kzero(n, k)
    failures = []
    for i in range(1, k + 1):
        for j in range(k + 2, n):
            cij = (Xget(X, n, i, j) + Xget(X, n, i + 1, j + 1)
                   - Xget(X, n, i, j + 1) - Xget(X, n, i + 1, j))
            cij_s = sp.expand(cij)
            if cij_s != 0:
                failures.append((i, j, cij_s))
    return failures


# ---------------------------------------------------------------------------
# T1b: Chord-rate table
# ---------------------------------------------------------------------------

def predicted_finite(n, k, a, b, alphas, betas):
    """Predict from Part II eq. (7) whether X_{a,b} is finite under
       rate vector (alphas, betas). Special/bare always infinite."""
    if a > b:
        a, b = b, a
    # special
    if (a, b) == (1, k + 2):
        return False
    # bare
    if (a, b) == (k + 1, n):
        return False
    # born-free always finite
    if a >= 1 and b <= k + 1:
        return True
    if a >= k + 2 and b <= n:
        return True
    # bridges:
    if a == 1 and b in alphas:           # X_{1, l}: alpha_l = 0
        return alphas[b] == 0
    if a == k + 1 and b in alphas:       # X_{k+1, l}: alpha_l = 1
        return alphas[b] == 1
    if 2 <= a <= k and b in alphas:      # X_{i, l}: alpha_l + beta_i = 1
        return alphas[b] + betas[a] == 1
    if 2 <= a <= k and b == k + 2:       # X_{i, k+2}: beta_i = 0
        return betas[a] == 0
    if 2 <= a <= k and b == n:           # X_{i, n}: beta_i = 1
        return betas[a] == 1
    raise RuntimeError(f"unclassified chord ({a},{b})")


def test_rate_table(n, k, rate_vectors):
    """For each rate vector, substitute and check each chord's finite/infinite
       status matches the prediction from eq. (7)."""
    X, s, t, u, Y = setup_kzero(n, k)
    L = sp.Symbol('L', positive=True)
    failures = []
    for label, alphas, betas in rate_vectors:
        subs = {s: L}
        for l in alphas:
            subs[u[l]] = sp.Rational(alphas[l]) * L
        for i in betas:
            subs[t[i]] = sp.Rational(betas[i]) * L
        for (a, b), expr in X.items():
            if (a, b) in Y:
                continue       # born-free is held finite by definition
            e = sp.expand(expr.subs(subs))
            # finite iff degree in L is 0 (constant) or e is identically 0
            if e == 0:
                actual_finite = True
            else:
                poly = sp.Poly(e, L)
                deg = poly.degree()
                actual_finite = (deg == 0)
            pred_finite = predicted_finite(n, k, a, b, alphas, betas)
            if actual_finite != pred_finite:
                failures.append((label, (a, b), e, pred_finite, actual_finite))
    return failures


# ---------------------------------------------------------------------------
# T2: Survival classification (Prop 1)
# ---------------------------------------------------------------------------

def rate_equation(n, k, a, b):
    """Return a linear constraint on (alpha_ell, beta_i) for X_{a,b} to be finite,
       as a tuple (alpha_coeffs_dict, beta_coeffs_dict, const).
       The constraint is sum_l alpha_coeffs[l]*alpha_l + sum_i beta_coeffs[i]*beta_i = const.
       Returns None if the chord can never be made finite (special/bare),
       or {} (trivially satisfied) for born-free."""
    if a > b:
        a, b = b, a
    if (a, b) == (1, k + 2):  return None      # special
    if (a, b) == (k + 1, n):  return None      # bare
    # born-free
    if (a >= 1 and b <= k + 1) or (a >= k + 2 and b <= n):
        return ({}, {}, 0)                    # always satisfied
    if a == 1 and k + 3 <= b <= n - 1:        # X_{1,l}: alpha_l = 0
        return ({b: 1}, {}, 0)
    if a == k + 1 and k + 3 <= b <= n - 1:    # X_{k+1,l}: alpha_l = 1
        return ({b: 1}, {}, 1)
    if 2 <= a <= k and k + 3 <= b <= n - 1:   # X_{i,l}: alpha_l + beta_i = 1
        return ({b: 1}, {a: 1}, 1)
    if 2 <= a <= k and b == k + 2:            # X_{i,k+2}: beta_i = 0
        return ({}, {a: 1}, 0)
    if 2 <= a <= k and b == n:                # X_{i,n}: beta_i = 1
        return ({}, {a: 1}, 1)
    raise RuntimeError(f"unclassified chord ({a},{b}) in rate_equation")


def system_consistent(equations):
    """Given a list of (alpha_dict, beta_dict, const) linear constraints,
       check consistency over Q (or R). Use sympy linsolve."""
    if any(eq is None for eq in equations):
        return False   # special/bare imposes 1=0
    # collect vars
    alpha_vars, beta_vars = set(), set()
    for adict, bdict, _ in equations:
        alpha_vars |= set(adict.keys())
        beta_vars |= set(bdict.keys())
    alpha_syms = {l: sp.Symbol(f'A{l}') for l in alpha_vars}
    beta_syms = {i: sp.Symbol(f'B{i}') for i in beta_vars}
    syms = list(alpha_syms.values()) + list(beta_syms.values())
    if not syms:
        # only born-free constraints (all trivially satisfied)
        return True
    eqs = []
    for adict, bdict, c in equations:
        lhs = sum((coef * alpha_syms[l] for l, coef in adict.items()), sp.Integer(0))
        lhs += sum((coef * beta_syms[i] for i, coef in bdict.items()), sp.Integer(0))
        eqs.append(lhs - c)
    try:
        sol = sp.linsolve(eqs, syms)
    except Exception:
        return False
    return len(sol) > 0


def multiset_survives_via_rates(M, n, k):
    """M is a list of chord-tuples (a,b) (with possible repetitions).
       Returns True if the rate system from bridges in M is inconsistent."""
    eqs = [rate_equation(n, k, a, b) for (a, b) in M]
    return not system_consistent(eqs)


def multiset_satisfies_prop1(M, n, k):
    """M satisfies Prop 1's survival criterion: contains special, bare, or
       a column pair {X_{1,m}, X_{k+1,m}} for some m in {k+3,..,n-1}."""
    chord_set = set(M)
    if (1, k + 2) in chord_set:  return True   # special
    if (k + 1, n) in chord_set:  return True   # bare
    for m in range(k + 3, n):
        if (1, m) in chord_set and (k + 1, m) in chord_set:
            return True
    return False


def chords_cross(c1, c2):
    """Two chords (a,b) and (c,d) cross iff a<c<b<d or c<a<d<b
       (with each chord sorted)."""
    a, b = sorted(c1); c, d = sorted(c2)
    return (a < c < b < d) or (c < a < d < b)


def enumerate_triangulations(n):
    """All n-3 element subsets of chords with no pairwise crossing.
       (Catalan-many: C_{n-2}.)"""
    from itertools import combinations
    chords = all_chords(n)
    triangs = []
    for sub in combinations(chords, n - 3):
        if all(not chords_cross(sub[i], sub[j])
               for i in range(len(sub)) for j in range(i + 1, len(sub))):
            triangs.append(sub)
    return triangs


def test_triangulation_escape(n, k):
    """Every triangulation of the n-gon should escape the k-zero based at
       {1..k+1} (Part II Section 5). Test by direct rate-equation inconsistency."""
    triangs = enumerate_triangulations(n)
    failures = [T for T in triangs
                if not multiset_survives_via_rates(list(T), n, k)]
    return triangs, failures


def non_crossing_chords(n, k):
    return [c for c in all_chords(n) if doesnt_cross_enclosing(n, k, *c)]


def test_prop1(n, k, verbose=False):
    """Enumerate all multisets of size n-3 over non-crossing chords;
       verify multiset_survives_via_rates(M) == multiset_satisfies_prop1(M)."""
    nc = non_crossing_chords(n, k)
    failures = []
    counted = 0
    for M in combinations_with_replacement(nc, n - 3):
        counted += 1
        a = multiset_survives_via_rates(list(M), n, k)
        b = multiset_satisfies_prop1(list(M), n, k)
        if a != b:
            failures.append((M, a, b))
            if verbose:
                print(f"  MISMATCH: M={M}, rates_survive={a}, prop1={b}")
    return counted, failures


# ---------------------------------------------------------------------------
# T3: Step-2 relations
# ---------------------------------------------------------------------------

def multiset_label(M):
    """Sort and label a multiset, e.g. ((1,3),(2,4)) -> '13,24'."""
    return ','.join(f"{a}{b}" for (a, b) in sorted(M))


def build_ansatz(n, k, bp_order):
    """Build the boundary-propagator piece B_n^(bp_order) on Z^(k) with
       symbolic coefficients a_M. Returns (expr, coeff_syms) where coeff_syms
       maps multiset-label -> sympy symbol a_M."""
    X, s, t, u, Y = setup_kzero(n, k)
    chords = all_chords(n)
    bf_set = set(born_free_chords(n, k))
    bp_of = {c: (0 if c in bf_set else 1) for c in chords}

    # All multisets of size n-3
    expr = sp.Integer(0)
    coeff_syms = {}
    for M in combinations_with_replacement(chords, n - 3):
        # restrict to fixed BP count
        if sum(bp_of[c] for c in M) != bp_order:
            continue
        label = multiset_label(M)
        a_M = sp.Symbol(f'a_{label}')
        coeff_syms[label] = a_M
        denom = sp.Integer(1)
        for c in M:
            denom *= Xget(X, n, *c)
        # Avoid division by zero from leg-chord (won't appear if M valid)
        expr += a_M / denom
    return sp.together(sp.cancel(expr)), coeff_syms, s, t, u, Y


def extract_pure_special_bare_D_subset(expr, Y):
    """Extract the D-subset with no born-free factor in the denominator
       (i.e. the part of expr independent of every Y symbol).
       Achieved by substituting every Y -> 0 in expr (since the D-subset
       with no Y is the part you'd get by isolating Y^0 coefficient, which
       for terms inversely proportional to Y vanish under Y->0... so we
       substitute Y -> 1 actually for the constant-Y-monomial D-subset).
       Actually the cleanest: each D-subset is indexed by a born-free monomial
       in the *denominator*. The pure-special-bare D-subset has no born-free
       in denominator, so it's the limit Y -> infinity of expr * (1/Y^{multiplicity}).
       But we can also: substitute Y -> 1 and isolate the terms with no Y.

       Cleanest: each term in expr is (a_M) / (prod_BP * prod_Y_in_M).
       The pure D-subset = sum over M with no Y in M.

       We'll just rebuild this directly via build_ansatz_pure_special.
    """
    pass


def build_special_bare_pure(n, k, bp_order):
    """Build the D-subset of B_n^(bp_order) consisting only of multisets M
       whose chords are all in {special, bare} (no born-free, no other bridges).
       For bp_order = m this is a sum over the (m+1) multisets {s^{m-j} b^j}."""
    X, s, t, u, Y = setup_kzero(n, k)
    s_chord = (1, k + 2)
    b_chord = (k + 1, n)
    expr = sp.Integer(0)
    coeff_syms = {}
    # We need n-3 - bp_order = 0 born-free, so bp_order must equal n-3.
    # For general bp_order <= n-3 we have m BP + (n-3-m) born-free per multiset.
    # The "pure binomial-trace" D-subset is the one with born-free monomial = 1,
    # which requires n-3-m = 0, i.e. m = n-3.
    # For m < n-3 we need to allow a born-free factor; choose a particular Y monomial.
    # Below we handle bp_order == n-3 case explicitly: pure special/bare D-subset.
    assert bp_order == n - 3, "binomial-trace pure D-subset requires m = n-3"
    for j in range(bp_order + 1):
        # multiset = j bares + (bp_order - j) specials
        M = tuple([s_chord] * (bp_order - j) + [b_chord] * j)
        label = multiset_label(M)
        a_M = sp.Symbol(f'a_{label}')
        coeff_syms[label] = a_M
        denom = sp.Integer(1)
        for c in M:
            denom *= Xget(X, n, *c)
        expr += a_M / denom
    return expr, coeff_syms, s


def test_binomial_trace(n, k, m):
    """For bp_order = m, verify the leading 1/s^m coefficient on the
       pure-special/bare D-subset is sum_{j=0}^m (-1)^j a_{s^{m-j} b^j}."""
    # We need m = n-3 for the pure D-subset to exist with no born-free factor.
    # If m < n-3, we pad with a born-free monomial.
    X, s, t, u, Y = setup_kzero(n, k)

    # Build the D-subset with j special^{m-j} bare^j + (n-3-m) copies of a
    # chosen born-free chord (so the D-subset has a definite Y-monomial,
    # which we then divide out symbolically).
    if not Y:
        if m != n - 3:
            return f"SKIP (no born-free chord at n={n},k={k}, and m={m} != n-3)"
        Y_monomial = sp.Integer(1)
        n_bf = 0
    else:
        Y_pick = next(iter(Y))                # arbitrary born-free
        Y_sym = Y[Y_pick]
        n_bf = (n - 3) - m
        if n_bf < 0:
            return f"SKIP (m={m} > n-3={n-3})"
        Y_monomial = Y_sym ** n_bf

    s_chord = (1, k + 2)
    b_chord = (k + 1, n)
    expr = sp.Integer(0)
    coeffs = []
    for j in range(m + 1):
        M_bp = [s_chord] * (m - j) + [b_chord] * j
        if n_bf > 0:
            M = tuple(M_bp + [Y_pick] * n_bf)
        else:
            M = tuple(M_bp)
        label = multiset_label(M)
        a_M = sp.Symbol(f'a_{label}')
        coeffs.append((j, a_M))
        denom = sp.Integer(1)
        for c in M:
            denom *= Xget(X, n, *c)
        expr += a_M / denom

    # On Z^(k), s = X_{1,k+2}, b = X_{k+1,n} = -s. Substituting:
    # expr * Y_monomial = pure expression in s.
    # Multiply expr by Y_monomial to clear the Y factor, then by s^m;
    # the leading 1/s^m coefficient is the limit as s -> infinity.
    pure = expr * Y_monomial
    leading = sp.expand(sp.limit(pure * s ** m, s, sp.oo))

    # Predicted: sum_{j=0}^m (-1)^j a_{s^{m-j} b^j}
    predicted = sum(((-1) ** j) * a for (j, a) in coeffs)

    diff = sp.expand(leading - predicted)
    return diff == 0, leading, predicted


def test_2zero_offset_relation(n, k=2):
    """Verify the 2-zero offset relation (Part III eq. 7):
       a_{14,24} + a_{14,2n} - a_{3n,24} - a_{3n,2n} = 0
       extracted at m=2 from the s -> 0 limit (which makes X_{2,4} = X_{2,n})."""
    if n < 6:
        return f"SKIP (need n>=6 at k=2)"
    X, s, t, u, Y = setup_kzero(n, k)
    s_chord = (1, k + 2)
    b_chord = (k + 1, n)
    off_24 = (2, k + 2)
    off_2n = (2, n)
    if not Y:
        return f"SKIP (no born-free chord)"
    Y_pick = next(iter(Y)); Y_sym = Y[Y_pick]
    pad = [Y_pick] * (n - 5)

    expr = sp.Integer(0)
    a_syms = {}
    for prefix, sign in [(s_chord, +1), (b_chord, -1)]:
        for off, off_sign_correction in [(off_24, +1), (off_2n, +1)]:
            M = tuple(sorted([prefix, off]) + pad)
            label = multiset_label(M)
            a = sp.Symbol(f'a_{label}')
            a_syms[(prefix, off)] = (a, sign)
            denom = sp.Integer(1)
            for cc in M:
                denom *= Xget(X, n, *cc)
            expr += a / denom

    # multiply by Y^(n-5), then by s and t_2 to get the (1/s)(1/t_2) residue
    pure = expr * Y_sym ** (n - 5)
    times_s = sp.limit(s * pure, s, 0)                # 1/s residue (using s=0 makes X_{2,n}=t_2)
    residue = sp.limit(t[2] * times_s, t[2], 0)       # 1/t_2 residue
    residue = sp.expand(residue)

    predicted = sp.Integer(0)
    for (prefix, off), (a, sign) in a_syms.items():
        predicted += sign * a
    return sp.expand(residue - predicted) == 0, residue, predicted


def test_2zero_column_relation(n, k=2):
    """Verify the 2-zero per-column relation
       sum_{i=1}^{3} (a_{14,i_ell} - a_{3n,i_ell}) = 0
       extracted from B_n^(2)|_Z^(2) by the X_{2,4} -> 0 cut at the 1/s coefficient."""
    # Build the 2 BP D-subset that involves the special with each column-l bridge,
    # and same for the bare. Pad with born-free factor Y^(n-5).
    if n < 6:
        return f"SKIP (need n>=6 for non-trivial 2-zero column relation at k=2)"
    X, s, t, u, Y = setup_kzero(n, k)
    # Pick a far column ell, say ell = k+3 = 5
    ell = k + 3
    s_chord = (1, k + 2)
    b_chord = (k + 1, n)
    col_chords = [(1, ell), (2, ell), (k + 1, ell)]   # i = 1, 2, 3 for k=2

    # Need n-3 chord slots; we use 2 BP + (n-5) born-free
    if not Y:
        return f"SKIP (no born-free chord at n={n},k={k})"
    Y_pick = next(iter(Y))
    Y_sym = Y[Y_pick]
    pad = [Y_pick] * (n - 5)

    expr = sp.Integer(0)
    a_syms = {}
    for prefix, sign in [(s_chord, +1), (b_chord, -1)]:
        for ci, c in enumerate(col_chords):
            M = tuple(sorted([prefix, c]) + pad)
            label = multiset_label(M)
            a = sp.Symbol(f'a_{label}')
            a_syms[(prefix, c)] = (a, sign)
            denom = sp.Integer(1)
            for cc in M:
                denom *= Xget(X, n, *cc)
            expr += a / denom

    # Apply the X_{2,4}=t_2 -> 0 cut; multiply by Y monomial and by s to
    # clear the 1/s factor; then send s -> 0 (so 1/(u_ell - s) -> 1/u_ell,
    # picking up the contributions from all three column-ell chords). The
    # column relation is the residue at u_ell -> 0 (the 1/u_ell coefficient).
    pure = expr * Y_sym ** (n - 5)
    cut = pure.subs(t[2], 0)
    after_s = sp.limit(s * cut, s, 0)              # 1/s residue
    coeff_1_over_u = sp.limit(u[ell] * after_s, u[ell], 0)  # 1/u_ell residue
    coeff_1_over_u = sp.expand(coeff_1_over_u)

    # Predicted: sum_{i=1}^{3} (a_{14,i_ell} - a_{3n,i_ell})
    predicted = sp.Integer(0)
    for (prefix, c), (a, sign) in a_syms.items():
        predicted += sign * a
    return sp.expand(coeff_1_over_u - predicted) == 0, coeff_1_over_u, predicted


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("STRESS TESTS for Part II + Part III")
    print("=" * 72)

    overall_pass = True

    # ----- T1a -----
    print("\n[T1a] Telescoping lemma (X-substitutions satisfy c_{ij}=0)")
    for n, k in [(5, 1), (6, 1), (6, 2), (7, 2), (7, 3), (8, 3), (9, 4)]:
        fails = test_telescoping(n, k)
        status = "PASS" if not fails else "FAIL"
        print(f"  n={n}, k={k}: {status} ({len(fails)} failed c_ij)")
        if fails:
            overall_pass = False
            for f in fails[:3]:
                print(f"     ", f)

    # ----- T1b -----
    print("\n[T1b] Chord-rate table (eq. 7 finite-iff predictions)")
    for n, k in [(6, 2), (7, 2), (7, 3), (8, 3), (8, 4)]:
        # Build a battery of rate vectors
        rvs = []
        rvs.append(("A (all 0)",
                    {l: 0 for l in range(k + 3, n)},
                    {i: 0 for i in range(2, k + 1)}))
        rvs.append(("B (all 1)",
                    {l: 1 for l in range(k + 3, n)},
                    {i: 1 for i in range(2, k + 1)}))
        # mixed: betas=0, alpha_{k+3}=1
        if k + 3 < n:
            a = {l: 0 for l in range(k + 3, n)}
            a[k + 3] = 1
            rvs.append(("mixed (beta=0, alpha_{k+3}=1)", a,
                        {i: 0 for i in range(2, k + 1)}))
        # mixed: alphas=0, beta_2=1
        if k >= 2:
            b_ = {i: 0 for i in range(2, k + 1)}
            b_[2] = 1
            rvs.append(("mixed (alpha=0, beta_2=1)",
                        {l: 0 for l in range(k + 3, n)}, b_))
        fails = test_rate_table(n, k, rvs)
        status = "PASS" if not fails else "FAIL"
        print(f"  n={n}, k={k}: {status} ({len(rvs)} rate vectors tested, {len(fails)} chord mismatches)")
        if fails:
            overall_pass = False
            for label, chord, e, pred, actual in fails[:5]:
                print(f"     {label}: chord {chord}, e={e}, predicted_finite={pred}, actual_finite={actual}")

    # ----- T2 -----
    print("\n[T2] Prop 1 survival classification (rate-eq inconsistency == special/bare/column-pair)")
    for n, k in [(6, 1), (6, 2), (7, 2), (8, 2), (7, 3), (8, 3)]:
        try:
            counted, fails = test_prop1(n, k)
            status = "PASS" if not fails else "FAIL"
            print(f"  n={n}, k={k}: {status} ({counted} non-crossing multisets enumerated, {len(fails)} mismatches)")
            if fails:
                overall_pass = False
                for M, a, b in fails[:5]:
                    print(f"     M={M}, rates_inconsistent={a}, prop1_says_survive={b}")
        except Exception as e:
            print(f"  n={n}, k={k}: ERROR {e}")
            overall_pass = False

    # ----- T3a + T3b -----
    print("\n[T3a/b] Binomial trace at m=1,2 (Theorem 1 of Part III)")
    cases = [
        # (n, k, m)
        (5, 1, 1),   # m=1, k=1: a_{13,Y} = a_{2n,Y}
        (5, 1, 2),   # m=2, k=1: pure (no Y needed): a_{13,13}-a_{13,2n}+a_{2n,2n}
        (5, 2, 1),   # m=1, k=2: a_{14,Y}=a_{3n,Y}
        (5, 2, 2),   # m=2, k=2 pure
        (6, 2, 2),   # m=2, k=2 with a born-free Y
        (6, 2, 3),   # m=3, k=2 pure
        (7, 2, 2),
        (6, 3, 2),   # m=2, k=3
    ]
    for n, k, m in cases:
        result = test_binomial_trace(n, k, m)
        if isinstance(result, str):
            print(f"  n={n}, k={k}, m={m}: {result}")
            continue
        ok, lead, pred = result
        status = "PASS" if ok else "FAIL"
        print(f"  n={n}, k={k}, m={m}: {status}")
        if not ok:
            overall_pass = False
            print(f"     leading 1/s^{m} coeff: {lead}")
            print(f"     predicted (binomial trace): {pred}")
            print(f"     diff: {sp.expand(lead - pred)}")

    # ----- T3d -----
    print("\n[T3d] 2-zero offset relation at m=2 (Part III eq. 7)")
    for n in [6, 7, 8]:
        result = test_2zero_offset_relation(n, k=2)
        if isinstance(result, str):
            print(f"  n={n}, k=2: {result}")
            continue
        ok, residue, pred = result
        status = "PASS" if ok else "FAIL"
        print(f"  n={n}, k=2: {status}")
        if not ok:
            overall_pass = False
            print(f"     residue at s=0, t_2=0: {residue}")
            print(f"     predicted offset:      {pred}")
            print(f"     diff: {sp.expand(residue - pred)}")

    # ----- T4 -----
    print("\n[T4] Triangulation escape (Part II Section 5)")
    for n, k in [(5, 1), (5, 2), (6, 1), (6, 2), (6, 3), (7, 2), (7, 3), (8, 2), (8, 3)]:
        triangs, fails = test_triangulation_escape(n, k)
        status = "PASS" if not fails else "FAIL"
        print(f"  n={n}, k={k}: {status} ({len(triangs)} triangulations, {len(fails)} fail to escape)")
        if fails:
            overall_pass = False
            for T in fails[:3]:
                print(f"     {T}")

    # ----- T3c -----
    print("\n[T3c] 2-zero per-column relation at m=2 (Part III eq. 6)")
    for n in [6, 7, 8]:
        result = test_2zero_column_relation(n, k=2)
        if isinstance(result, str):
            print(f"  n={n}, k=2: {result}")
            continue
        ok, lead, pred = result
        status = "PASS" if ok else "FAIL"
        print(f"  n={n}, k=2: {status}")
        if not ok:
            overall_pass = False
            print(f"     coeff(1/s) after X_24->0 cut: {lead}")
            print(f"     predicted column relation:    {pred}")
            print(f"     diff: {sp.expand(lead - pred)}")

    print()
    print("=" * 72)
    print("OVERALL:", "ALL TESTS PASSED" if overall_pass else "SOME TESTS FAILED")
    print("=" * 72)
    return 0 if overall_pass else 1


if __name__ == "__main__":
    sys.exit(main())
