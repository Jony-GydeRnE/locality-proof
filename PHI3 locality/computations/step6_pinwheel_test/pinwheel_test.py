#!/usr/bin/env python3
"""
step6_pinwheel_test/pinwheel_test.py
================================================================================

Tests the BOUNDED-K LOCALITY conjecture by building m-fold symmetric pinwheels
at n = 3m and checking their minimal killing k under the GEOMETRIC P1/P2 leading
kill (Part II Theorem 1).

Conjecture (Jonathan + first-Claude analysis, 2026-05-31):
    K(n) := max over non-tri multisets M at the n-gon of (min k that kills M
            under the geometric leading kill at some rotation i)
    grows like K(n) = floor(n/3).
    Data points: K(7) = K(8) = 2, K(9) = 3.

The 57 multisets at n=9 that need k=3 (vs k=2 for the other 954 step-1
survivors) all contain a TIGHT crossing pair -- two short chords (j,j+2) and
(j+1,j+3) crossing on 4 consecutive vertices. The most symmetric of these is
the Z_3 pinwheel at n=9: three tight crossings rotated by n/3.

This script builds the analogous m-fold pinwheel at n = 3m for m=3,4,5 and
checks: does the m-fold pinwheel actually need k=m to be killed?

KILL CRITERION
--------------
GEOMETRIC (matches step1_uncaught.py at k=1; matches first-Claude's
kill_dictionary.py framework):
    M survives the k-zero at rotation i iff:
        (P1) some chord of M crosses the enclosing chord X_{i, i+k}, OR
        (P2) M contains the special X_{i, i+k+1}, the bare X_{i+k, i-1} (cyclic),
             or a column pair {X_{i, m}, X_{i+k, m}} for some far middle m.
    Killed at (k, i) iff neither P1 nor P2 holds.
    At k=1, P1 is vacuous (no real enclosing chord), so this reduces to the
    step-1 criterion (no special / no bare / no column pair in M).

NB: this geometric criterion is OVER-PERMISSIVE relative to the rigorous
rate-equation kill (a single P1 crossing chord can typically be held finite
under a careful rate choice). Under the rigorous criterion every step-1
survivor at n=9 dies at k=2 (see step5 findings). The geometric K(n) >= 3 at
n=9 is a feature of the geometric heuristic and reveals real combinatorial
structure (the tight-crossing pinwheel) even though it disagrees with the
algebraic kill.

USAGE
-----
    python3 pinwheel_test.py              # validate at n=9 + run pinwheel m=3,4,5
    python3 pinwheel_test.py --only-pinwheel  # skip n=9 validation
"""
import argparse
import json
import os
import sys
from itertools import combinations_with_replacement, combinations


# =========================================================================
# 1. Chord / multiset machinery
# =========================================================================

def normalize(i, j, n):
    """Canonical chord (a,b) with a<b, a,b in {1..n}, after cyclic reduction."""
    i = ((i - 1) % n) + 1
    j = ((j - 1) % n) + 1
    if i > j:
        i, j = j, i
    return (i, j)


def is_leg(n, a, b):
    if a > b:
        a, b = b, a
    return (b - a == 1) or (a == 1 and b == n)


def all_chords(n):
    return [(i, j) for i in range(1, n + 1)
                   for j in range(i + 2, n + 1)
                   if not (i == 1 and j == n)]


def chords_cross(c1, c2):
    """Two chords cross iff their endpoints alternate around the boundary."""
    a, b = c1
    c, d = c2
    if {a, b} & {c, d}:
        return False
    inside = lambda v: a < v < b
    return inside(c) != inside(d)


def is_triangulation(ms, n):
    if len(set(ms)) != len(ms):
        return False
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            if chords_cross(ms[i], ms[j]):
                return False
    return True


def chord_crosses_chord(chord, enclosing_chord, n):
    """Does `chord` cross `enclosing_chord` in the n-gon?
    Treat enclosing as a real chord (k>=2). For k=1 the 'enclosing chord' is an
    edge and is degenerate; in that case nothing crosses it geometrically."""
    if enclosing_chord is None:
        return False
    return chords_cross(chord, enclosing_chord)


# =========================================================================
# 2. Geometric kill criterion at (k, i)
# =========================================================================

def kzone_data(n, k, i):
    """For the k-zero based at vertices {i, i+1, ..., i+k}:
        enclosing : chord X_{i, i+k}; None at k=1 (it's an edge).
        special   : X_{i, i+k+1}.
        bare      : X_{i+k, i-1} cyclic.
        column_pairs : list of (X_{i, m}, X_{i+k, m}) for far-middle m.
    """
    if k == 1:
        enclosing = None        # X_{i, i+1} is an edge
    else:
        enclosing = normalize(i, i + k, n)
    special = normalize(i,     i + k + 1, n)
    bare    = normalize(i + k, i - 1,     n)
    # Far middle vertices: {i+k+1, ..., i-1} cyclic excluding the special-far
    # neighbor i+k+1 itself and the bare-far neighbor i-1.
    near = set(((i - 1 + a) % n) + 1 for a in range(k + 1))
    far  = [v for v in range(1, n + 1) if v not in near]
    sp_far = ((i + k + 1 - 1) % n) + 1
    br_far = ((i - 1 - 1) % n) + 1
    column_pairs = []
    for m in far:
        if m == sp_far or m == br_far:
            continue
        c_low  = normalize(i,     m, n)
        c_high = normalize(i + k, m, n)
        if c_low == c_high:
            continue
        if is_leg(n, *c_low) or is_leg(n, *c_high):
            continue
        column_pairs.append((c_low, c_high))
    return {'k': k, 'i': i, 'n': n,
            'enclosing': enclosing,
            'special': special, 'bare': bare,
            'column_pairs': column_pairs}


def survives_geometric_at(ms_set, kz):
    """M survives the (k, i) zone under the geometric criterion."""
    # P2 first (cheap)
    if kz['special'] in ms_set:  return True
    if kz['bare']    in ms_set:  return True
    for (c1, c2) in kz['column_pairs']:
        if c1 in ms_set and c2 in ms_set:
            return True
    # P1 -- any chord of M crosses the enclosing chord
    enc = kz['enclosing']
    if enc is not None:
        for c in ms_set:
            if chords_cross(c, enc):
                return True
    return False


def killed_by_some_zone_geom(ms, n):
    """Return (True, (k, i)) for the lex-min (k, i) that geometrically kills M,
    or (False, None) if M survives every (k, i)."""
    ms_set = set(ms)
    for k in range(1, n - 2):
        for i in range(1, n + 1):
            kz = kzone_data(n, k, i)
            if not survives_geometric_at(ms_set, kz):
                return True, (k, i)
    return False, None


def min_killing_k_geom(ms, n):
    """Return the minimum k such that M is killed at (k, i) for some i, under
    the geometric criterion. Returns None if M survives every (k, i)."""
    ms_set = set(ms)
    for k in range(1, n - 2):
        for i in range(1, n + 1):
            kz = kzone_data(n, k, i)
            if not survives_geometric_at(ms_set, kz):
                return k, i
    return None, None


# =========================================================================
# 3. Validate at n=9: confirm the geometric criterion gives the same
#    954@k=2, 57@k=3 split as the first-Claude analysis.
# =========================================================================

def validate_n9():
    """Reproduce the first-Claude split: 1011 step-1 non-tri survivors,
    of which 954 die at k=2 (geom) and 57 need k=3 (geom)."""
    n = 9
    print(f"\n[validation at n={n}]  enumerating step-1 survivors and checking "
          f"min-kill-k under the GEOMETRIC criterion...")
    chords = all_chords(n)
    N = n - 3
    step1_survivors = []
    for ms in combinations_with_replacement(chords, N):
        if is_triangulation(ms, n):
            continue
        killed_at_k1 = False
        for i in range(1, n + 1):
            kz = kzone_data(n, 1, i)
            if not survives_geometric_at(set(ms), kz):
                killed_at_k1 = True
                break
        if not killed_at_k1:
            step1_survivors.append(ms)
    n_total = len(step1_survivors)

    # For each, find min-kill-k under geometric criterion across k=2..n-3
    by_k = {}
    examples_by_k = {}
    for ms in step1_survivors:
        ms_set = set(ms)
        min_k = None
        min_ki = None
        for k in range(2, n - 2):
            for i in range(1, n + 1):
                kz = kzone_data(n, k, i)
                if not survives_geometric_at(ms_set, kz):
                    min_k = k
                    min_ki = (k, i)
                    break
            if min_k is not None:
                break
        by_k[min_k] = by_k.get(min_k, 0) + 1
        if min_k not in examples_by_k:
            examples_by_k[min_k] = []
        if len(examples_by_k.get(min_k, [])) < 3:
            examples_by_k[min_k].append((ms, min_ki))

    print(f"  step-1 non-tri survivors at n={n}: {n_total}")
    print(f"  min-killing-k histogram (under GEOMETRIC criterion):")
    for k, c in sorted(by_k.items(), key=lambda x: (x[0] is None, x[0])):
        print(f"    k = {k}: {c}")
    return by_k, examples_by_k


# =========================================================================
# 4. m-fold pinwheel construction at n = 3m
# =========================================================================

def pinwheel(n, m):
    """Build the m-fold pinwheel at n = 3m vertices.

    Each tight crossing uses 2 chords on 4 consecutive vertices
    {3j+1, 3j+2, 3j+3, 3j+4} cyclically. The crossing pair is the two
    short diagonals (3j+1, 3j+3) and (3j+2, 3j+4) (cyclic).

    Returns: list of 2m chord-pairs. Multiset size = 2m.
    """
    assert n == 3 * m, f"pinwheel requires n = 3m; got n={n}, m={m}"
    chords = []
    for j in range(m):
        v1 = ((3 * j + 0) % n) + 1
        v2 = ((3 * j + 1) % n) + 1
        v3 = ((3 * j + 2) % n) + 1
        v4 = ((3 * j + 3) % n) + 1
        chords.append(normalize(v1, v3, n))
        chords.append(normalize(v2, v4, n))
    return chords


def pinwheel_at_n9():
    """The Z_3 pinwheel at n=9 from first-Claude analysis."""
    return pinwheel(9, 3)


# =========================================================================
# 5. Test: build m-fold pinwheel + fillers, compute min-kill-k
# =========================================================================

def all_filler_choices(n, pinwheel_chords, num_fillers):
    """Return all possible filler-multisets (size num_fillers) where each
    chord is from the n-gon's chord set (born-free or pinwheel-overlapping
    chords allowed). Order-independent: combinations_with_replacement."""
    chords = all_chords(n)
    return combinations_with_replacement(chords, num_fillers)


def pinwheel_test(n, m, max_fillers_to_try=None, verbose=True):
    """Build the m-fold pinwheel at n = 3m. Need n-3 chords total; the pinwheel
    has 2m, so n-3 - 2m fillers are needed.

    For each choice of fillers, compute the min-kill-k.

    Returns: dict with histogram of min-kill-k across fillers, plus extremes.
    """
    pw = pinwheel(n, m)
    n_pw = len(pw)
    n_fillers = (n - 3) - n_pw
    assert n_fillers >= 0, f"pinwheel has more chords than fits in size-{n-3} multiset"

    if verbose:
        print(f"\n[pinwheel n={n}, m={m}]")
        print(f"  pinwheel chords ({n_pw}): {pw}")
        print(f"  fillers needed: {n_fillers}")

    by_k = {}
    examples_by_k = {}
    n_filler_choices = 0
    fillers_iter = all_filler_choices(n, pw, n_fillers)
    for filler in fillers_iter:
        # full multiset = pinwheel + filler, sorted
        M = tuple(sorted(list(pw) + list(filler)))
        n_filler_choices += 1
        if max_fillers_to_try is not None and n_filler_choices > max_fillers_to_try:
            break
        # Some fillers may make M a triangulation -- skip those (triangulations
        # always survive every (k,i); we want the min-kill-k for the survivors
        # but a tri 'survives' all zones via P2 cleanly).
        # NOTE: tri's "survive" geometrically because some completion always
        # holds. We're interested in the non-tri case primarily.
        is_tri = is_triangulation(M, n)
        min_k_pair = min_killing_k_geom(M, n)
        min_k = min_k_pair[0]      # k or None
        key = (min_k, is_tri)
        by_k[key] = by_k.get(key, 0) + 1
        if key not in examples_by_k:
            examples_by_k[key] = []
        if len(examples_by_k[key]) < 3:
            examples_by_k[key].append({'M': M, 'min_killing_ki': min_k_pair})

    if verbose:
        print(f"  tested {n_filler_choices} filler choices")
        print(f"  min-kill-k histogram (key=(k_or_None, is_triangulation)):")
        for key, c in sorted(by_k.items(),
                             key=lambda x: ((x[0][0] is None), x[0][0] or 0,
                                            x[0][1])):
            k, tri = key
            label = f"k={k}" if k is not None else "SURVIVES ALL"
            tri_label = " (triangulation)" if tri else ""
            print(f"    {label}{tri_label}: {c}")

    return {
        'n': n, 'm': m,
        'pinwheel_chords': pw,
        'num_fillers': n_fillers,
        'n_filler_choices_tested': n_filler_choices,
        'min_k_histogram': by_k,
        'examples_by_k': examples_by_k,
    }


# =========================================================================
# 6. Pure-pinwheel test: when the pinwheel exactly fills the multiset slots
#    (n = 3m with n-3 = 2m, i.e., m = 3), no fillers needed. Test directly.
# =========================================================================

def pure_pinwheel_min_k(n, m):
    """For the exact-fit case (n = 3m, multiset size 2m = n-3 means n = 3 ???).
    Actually n-3 = 2m only at m=3 (n=9): 2(3)=6, n-3=6. So at n=9 m=3 we have
    exact fit; for m>=4 we need fillers."""
    pw = pinwheel(n, m)
    if len(pw) != n - 3:
        return None    # not exact fit
    M = tuple(sorted(pw))
    return min_killing_k_geom(M, n)


# =========================================================================
# 7. Main
# =========================================================================

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only-pinwheel', action='store_true',
                    help='Skip n=9 validation, jump straight to pinwheel test')
    ap.add_argument('--cap-fillers', type=int, default=None,
                    help='Cap the number of filler choices tested at each n')
    ap.add_argument('--max-n', type=int, default=15,
                    help='Largest n = 3m to test (default 15 = m=5)')
    args = ap.parse_args()

    results = {}

    # ----- Validation at n = 9 -----
    if not args.only_pinwheel:
        by_k, examples_by_k = validate_n9()
        results['validation_n9'] = {
            'min_k_histogram': by_k,
            'expected': '954 at k=2, 57 at k=3 per first-Claude',
        }
        # Save
        with open(os.path.join(OUT_DIR, 'validation_n9.json'), 'w') as f:
            json.dump({
                'min_k_histogram': {str(k): v for k, v in by_k.items()},
            }, f, indent=2)

    # ----- Exact-fit pinwheel at n=9 (m=3) -----
    print(f"\n[exact pinwheel] n=9, m=3 -- the Z_3 pinwheel with no fillers")
    M_pw9 = tuple(sorted(pinwheel(9, 3)))
    print(f"  M = {M_pw9}")
    k_kill, ki = min_killing_k_geom(M_pw9, 9)
    print(f"  min killing k = {k_kill} (at zone {ki})")
    results['exact_pinwheel_n9'] = {'M': M_pw9, 'min_k': k_kill, 'min_ki': ki}

    # ----- Pinwheel test at n=12, m=4 and n=15, m=5 -----
    for m in range(4, args.max_n // 3 + 1):
        n = 3 * m
        if n > args.max_n:
            break
        res = pinwheel_test(n, m, max_fillers_to_try=args.cap_fillers)
        results[f'pinwheel_n{n}_m{m}'] = {
            'pinwheel_chords': res['pinwheel_chords'],
            'num_fillers': res['num_fillers'],
            'n_filler_choices_tested': res['n_filler_choices_tested'],
            'min_k_histogram': {
                str(k): v for k, v in res['min_k_histogram'].items()
            },
        }
        # Pull out the worst-case (max min-kill-k) example
        non_tri_keys = [k for k in res['min_k_histogram'].keys() if not k[1]]
        non_tri_ks = [k[0] for k in non_tri_keys if k[0] is not None]
        if non_tri_ks:
            max_min_k = max(non_tri_ks)
            print(f"  WORST-CASE (max min-kill-k among non-tri): k = {max_min_k}")
            results[f'pinwheel_n{n}_m{m}']['max_min_k_non_tri'] = max_min_k
        if any(k[0] is None and not k[1] for k in res['min_k_histogram'].keys()):
            n_survives = res['min_k_histogram'][(None, False)]
            print(f"  WARNING: {n_survives} non-tri multisets survive EVERY (k,i) "
                  f"under geometric criterion")

        # Save
        with open(os.path.join(OUT_DIR, f'pinwheel_n{n}_m{m}.json'), 'w') as f:
            json.dump({
                'n': n, 'm': m,
                'num_fillers': res['num_fillers'],
                'n_filler_choices_tested': res['n_filler_choices_tested'],
                'min_k_histogram': [
                    {'min_k': k[0], 'is_triangulation': k[1], 'count': v}
                    for k, v in res['min_k_histogram'].items()
                ],
                'examples_by_k': [
                    {
                        'min_k': key[0], 'is_triangulation': key[1],
                        'samples': [
                            {'M': [list(c) for c in e['M']],
                             'min_killing_ki': list(e['min_killing_ki'])
                                                if e['min_killing_ki'][0] is not None else None}
                            for e in ex
                        ],
                    }
                    for key, ex in res['examples_by_k'].items()
                ],
            }, f, indent=2)

    # ----- Save overall summary -----
    with open(os.path.join(OUT_DIR, 'summary.json'), 'w') as f:
        json.dump({k: (str(v) if isinstance(v, dict) and any(
                isinstance(kk, tuple) for kk in v.keys()) else v)
                  for k, v in results.items()}, f, indent=2, default=str)
    print(f"\nwrote summary.json + per-test JSONs into {OUT_DIR}/")


if __name__ == "__main__":
    main()
