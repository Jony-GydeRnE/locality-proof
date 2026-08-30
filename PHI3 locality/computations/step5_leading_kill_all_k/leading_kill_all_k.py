#!/usr/bin/env python3
"""
step5_leading_kill_all_k.py
================================================================================

DECISIVE TEST: does the leading fat-zero kill (Part II's survival criterion,
applied at every k and every cyclic rotation) ALONE force locality?

If yes (World A), locality is pure Part II combinatorial geometry; no rank
computation, no step-2, no induction is needed for locality. Step-2 (Part III)
is then only for unitarity.

If no (World B), the multisets that survive every leading kill are exactly
what step-2 must kill. Their count + structure tells us how much step-2 work
is actually needed.

Hypothesis (Jonathan, 2026-05-31): World A. The k=1 leading-kill survivors at
n=7,8,9 (the "fish" / 100 / 1011) die under k>=2 zones, and there are no
crossing-only survivors that escape every (k,i).

METHOD
------
For each (k, i) with k in {1..n-3} and i in {1..n}:
    The k-zero based at vertices {i, i+1, ..., i+k} imposes a linear
    constraint surface on the X_ij's. Pick a kill limit: the special chord
    X_{i,i+k+1} -> infinity. Other chords each have a forced finite-or-not
    status depending on which "rates" we pick for the other free coordinates.

    A multiset M is KILLABLE at (k, i) iff there exists a rate assignment
    making every chord in M finite. Equivalently, the rate-equation system
    (one linear constraint per chord in M) is feasible.

    A multiset M SURVIVES the (k, i) leading kill iff no rate assignment
    works -- the rate-equation system is infeasible.

PRIMARY CHECK: rate-equation feasibility (exact, via sympy linsolve).
CROSS-CHECK: geometric V1-V4 (special in M / bare in M / column pair / offset
pair). V1-V4 are SUFFICIENT for survival but not necessary in general; they
miss triple-and-higher chord conflicts that arise when M contains crossings
of the enclosing chord plus offsets.

We report disagreements between primary and cross-check as a sanity flag.

OUTPUT
------
outputs/
    n7_findings.md, n8_findings.md, n9_findings.md
    n7_phase0.json, n7_phase1.json (per-survivor minimal (k,i) kill data)
    ... same for n=8, n=9
    summary.md  -- the World A/B verdict
"""

import argparse
import json
import os
import sys
import time
from itertools import combinations, combinations_with_replacement
from fractions import Fraction
import sympy as sp


# =========================================================================
# 1. Chord / multiset machinery (copied + extended from step1_uncaught.py)
# =========================================================================

def normalize(i, j, n):
    """Canonical chord name (a,b) with a<b, both in {1..n}, after cyclic
    reduction. Allows negative indices and overflow."""
    i = ((i - 1) % n) + 1
    j = ((j - 1) % n) + 1
    if i > j:
        i, j = j, i
    return (i, j)


def is_leg(n, a, b):
    """X_{a,b} is a leg iff |b-a|=1 or {a,b}={1,n}."""
    if a > b:
        a, b = b, a
    return (b - a == 1) or (a == 1 and b == n)


def all_chords(n):
    """All planar chords -- pairs (i,j) with 1<=i<j<=n that are not edges."""
    return [(i, j) for i in range(1, n + 1)
                   for j in range(i + 2, n + 1)
                   if not (i == 1 and j == n)]


def chords_cross(c1, c2):
    """Two chords cross iff their endpoints alternate on the boundary."""
    a, b = c1
    c, d = c2
    if {a, b} & {c, d}:
        return False
    inside = lambda v: a < v < b
    return inside(c) != inside(d)


def is_triangulation(ms, n):
    """A multiset is a triangulation iff all distinct AND no pairwise crossing."""
    if len(set(ms)) != len(ms):
        return False
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            if chords_cross(ms[i], ms[j]):
                return False
    return True


# =========================================================================
# 2. Zone structure for any (k, i)
# =========================================================================

def kzone_structure(n, k, i):
    """For the k-zero based at vertices {i, i+1, ..., i+k} (cyclic, 1-indexed),
    return the structure needed for the geometric V1-V4 survival check.

    Returns dict with keys:
        special       : X_{i, i+k+1}                       (always infinite)
        bare          : X_{i+k, i-1} cyclic                (always infinite)
        column_pairs  : list of (X_{i, m}, X_{i+k, m}) for far-middle m
        offset_pairs  : list of (X_{p, i+k+1}, X_{p, i-1}) for interior near p
        near_block    : set of vertices {i, i+1, ..., i+k} (cyclic)
        far_block     : set of vertices {i+k+1, ..., i-1} (cyclic)
        enclosing     : the chord X_{i, i+k}
    """
    near = set(((i - 1 + a) % n) + 1 for a in range(k + 1))   # {i, i+1, ..., i+k}
    far  = set(v for v in range(1, n + 1) if v not in near)

    special   = normalize(i, i + k + 1, n)
    bare      = normalize(i + k, i - 1, n)
    enclosing = normalize(i, i + k, n) if not is_leg(n, i, i + k) else None

    # Column pairs: for each "middle" far vertex m, pair (X_{i,m}, X_{i+k,m})
    # The "middle" far vertices are those that are not the special-side neighbor
    # (i+k+1) and not the bare-side neighbor (i-1).
    column_pairs = []
    sp_far = ((i + k + 1 - 1) % n) + 1
    br_far = ((i - 1 - 1) % n) + 1
    for m in far:
        if m == sp_far or m == br_far:
            continue
        c_low  = normalize(i,     m, n)
        c_high = normalize(i + k, m, n)
        if c_low == c_high:
            continue          # shouldn't happen given m is far-middle
        if is_leg(n, *c_low) or is_leg(n, *c_high):
            continue
        column_pairs.append((c_low, c_high))

    # Offset pairs: for each interior near vertex p in {i+1, ..., i+k-1},
    # pair (X_{p, i+k+1}, X_{p, i-1})
    offset_pairs = []
    for p_off in range(1, k):
        p = ((i - 1 + p_off) % n) + 1
        c_left  = normalize(p, i + k + 1, n)
        c_right = normalize(p, i - 1, n)
        if is_leg(n, *c_left) or is_leg(n, *c_right):
            continue
        if c_left == c_right:
            continue
        offset_pairs.append((c_left, c_right))

    return {
        'k': k, 'i': i, 'n': n,
        'special': special,
        'bare': bare,
        'enclosing': enclosing,
        'near_block': frozenset(near),
        'far_block': frozenset(far),
        'column_pairs': column_pairs,
        'offset_pairs': offset_pairs,
    }


def survives_geometric(ms_set, zs):
    """V1-V4 survival check (sufficient but not necessary).

    Returns True if any of V1 (special in M), V2 (bare in M),
    V3 (column pair fully in M), V4 (offset pair fully in M) holds.
    """
    if zs['special'] in ms_set:
        return True
    if zs['bare'] in ms_set:
        return True
    for (c1, c2) in zs['column_pairs']:
        if c1 in ms_set and c2 in ms_set:
            return True
    for (c1, c2) in zs['offset_pairs']:
        if c1 in ms_set and c2 in ms_set:
            return True
    return False


# =========================================================================
# 3. Rate-equation feasibility (RIGOROUS criterion)
#
# For the canonical k-zero at i=1, the chord values on Z^(k) are:
#   X_{1,k+2}   = s           (special; always infinite)
#   X_{i,k+2}   = t_i          (offsets, i=2..k)
#   X_{1,l}     = u_l          (l=k+3..n-1)
#   X_{k+1,l}   = u_l - s
#   X_{i,l}     = u_l - s + t_i  (i=2..k, l=k+3..n-1)
#   X_{k+1,n}   = -s           (bare; always infinite)
#   X_{i,n}     = -s + t_i     (i=2..k)
#   born-free chords          = independent free coords, always finite
#
# A chord stays finite under the limit s -> infty iff a particular linear
# combination of (u_l, t_i, s) is bounded -- equivalently the "rate" of
# that combination equals 0. With rates a_l = rate(u_l), b_i = rate(t_i),
# we get one linear "finite-iff" constraint per non-born-free chord.
#
# A multiset M is KILLABLE at (k, i=1) iff the system of constraints from
# its chords is feasible. We use sympy's linsolve over Q.
# =========================================================================

def rate_constraint(n, k, a, b):
    """Linear constraint on (a_l, b_i) for X_{a,b} to be finite under the
    canonical k-zero limit. Returns:
        ({alpha_dict}, {beta_dict}, const)  -- meaning sum a_dict[l] * alpha_l
                                               + sum b_dict[i] * beta_i = const
        None                                -- chord can never be finite
        ({}, {}, 0)                         -- trivially satisfied (born-free)
    """
    if a > b:
        a, b = b, a
    if (a, b) == (1, k + 2):    return None      # special
    if (a, b) == (k + 1, n):    return None      # bare
    # Born-free: chords internal to N = {1..k+1} or to F = {k+2..n}
    if (1 <= a and b <= k + 1) or (k + 2 <= a and b <= n):
        return ({}, {}, 0)
    if a == 1 and (k + 3) <= b <= (n - 1):
        return ({b: 1}, {}, 0)                   # X_{1,l}: alpha_l = 0
    if a == k + 1 and (k + 3) <= b <= (n - 1):
        return ({b: 1}, {}, 1)                   # X_{k+1,l}: alpha_l = 1
    if 2 <= a <= k and (k + 3) <= b <= (n - 1):
        return ({b: 1}, {a: 1}, 1)               # X_{i,l}: alpha_l + beta_i = 1
    if 2 <= a <= k and b == k + 2:
        return ({}, {a: 1}, 0)                   # X_{i,k+2}: beta_i = 0
    if 2 <= a <= k and b == n:
        return ({}, {a: 1}, 1)                   # X_{i,n}: beta_i = 1
    raise RuntimeError(f"rate_constraint: unclassified chord ({a},{b}) at n={n}, k={k}")


def killable_rate_system(ms_canonical, n, k):
    """ms_canonical is the multiset M in canonical (i=1) form. Test if the
    rate-equation system imposed by M's chords is feasible (consistent). If
    feasible, M is killable at this (k, i=1) zone. If infeasible, M survives."""
    eqs = [rate_constraint(n, k, a, b) for (a, b) in ms_canonical]
    if any(e is None for e in eqs):
        return False                              # contains special or bare
    # Collect variables
    alpha_vars, beta_vars = set(), set()
    for adict, bdict, _ in eqs:
        alpha_vars |= set(adict.keys())
        beta_vars  |= set(bdict.keys())
    if not alpha_vars and not beta_vars:
        # No rate constraints -> trivially feasible (all chords born-free)
        return True
    alpha_syms = {l: sp.Symbol(f"A{l}") for l in alpha_vars}
    beta_syms  = {i: sp.Symbol(f"B{i}") for i in beta_vars}
    syms = list(alpha_syms.values()) + list(beta_syms.values())
    equations = []
    for adict, bdict, c in eqs:
        if not adict and not bdict and c == 0:
            continue                              # tautology, skip
        lhs = sum((coef * alpha_syms[l] for l, coef in adict.items()), sp.Integer(0))
        lhs += sum((coef * beta_syms[i]  for i, coef in bdict.items()), sp.Integer(0))
        equations.append(lhs - c)
    if not equations:
        return True
    try:
        sol = sp.linsolve(equations, syms)
    except Exception:
        return False
    return len(sol) > 0


def canonicalize_multiset(ms, n, i):
    """Cyclically shift the multiset M so that vertex i becomes vertex 1.
    Shift: v -> ((v - i) mod n) + 1. Returns the relabeled multiset."""
    out = []
    for (a, b) in ms:
        a2 = ((a - i) % n) + 1
        b2 = ((b - i) % n) + 1
        if a2 > b2:
            a2, b2 = b2, a2
        out.append((a2, b2))
    return out


def killable_at_kzone(ms, n, k, i):
    """Is multiset M killable at the k-zero based at vertex i?
    PRIMARY criterion: rate-equation feasibility (canonicalize, then test)."""
    ms_canon = canonicalize_multiset(ms, n, i)
    return killable_rate_system(ms_canon, n, k)


def killed_by_some_zone(ms, n):
    """M is killed by SOME (k, i) leading kill iff killable at some (k, i).
    Returns (True, (k, i)) for the first killing zone, or (False, None) if M
    survives every leading kill."""
    for k in range(1, n - 2):
        for i in range(1, n + 1):
            if killable_at_kzone(ms, n, k, i):
                return True, (k, i)
    return False, None


def killing_zones(ms, n):
    """Return the full list of (k, i) zones where M is killable."""
    out = []
    for k in range(1, n - 2):
        for i in range(1, n + 1):
            if killable_at_kzone(ms, n, k, i):
                out.append((k, i))
    return out


# =========================================================================
# 4. Cross-check: geometric V1-V4 vs rate-equation
# =========================================================================

def cross_check(ms, n, verbose=False):
    """For each (k, i), compare the geometric survives_geometric() with
    rate-equation killable_rate_system(). Any disagreement (in particular:
    rate says 'survives' but geometric says 'killed') is a sanity flag.
    Returns list of (k, i, rate_killed, geom_killed)."""
    ms_set = set(ms)
    discrepancies = []
    for k in range(1, n - 2):
        for i in range(1, n + 1):
            zs = kzone_structure(n, k, i)
            geom_surv = survives_geometric(ms_set, zs)
            rate_kill = killable_at_kzone(ms, n, k, i)
            geom_kill = not geom_surv
            if rate_kill != geom_kill:
                discrepancies.append((k, i, rate_kill, geom_kill))
                if verbose:
                    print(f"   DISCREPANCY at (k={k}, i={i}): rate killable={rate_kill}, "
                          f"geom killable={geom_kill}")
    return discrepancies


# =========================================================================
# 5. Self-check: reproduce step-1 survivor counts at n=7, 8
# =========================================================================

def selfcheck_step1(n, expected):
    """Verify: restricted to k=1 only, our code produces `expected` survivors."""
    print(f"\n[self-check k=1 only at n={n}]  expected step-1 survivors: {expected}")
    survivors = []
    for ms in combinations_with_replacement(all_chords(n), n - 3):
        killable_k1 = False
        for i in range(1, n + 1):
            if killable_at_kzone(ms, n, 1, i):
                killable_k1 = True
                break
        if not killable_k1:
            survivors.append(ms)
    n_tri    = sum(1 for ms in survivors if     is_triangulation(ms, n))
    n_nontri = sum(1 for ms in survivors if not is_triangulation(ms, n))
    print(f"  step-1 survivors at n={n}: {len(survivors)} total ({n_tri} tri, {n_nontri} non-tri)")
    assert n_nontri == expected, f"FAIL: expected {expected} non-tri step-1 survivors, got {n_nontri}"
    print(f"  -> matches expected {expected} non-tri survivors. OK.")
    return survivors


# =========================================================================
# 6. Phase 0 / Phase 1 runners
# =========================================================================

def phase0_for_n(n, step1_survivors_nontri):
    """Take known non-triangulation step-1 (k=1) survivors and check whether
    each survives every k>=2 zone (and also k=1, just to be sure)."""
    print(f"\n[Phase 0 at n={n}]  testing {len(step1_survivors_nontri)} step-1 non-tri survivors "
          f"against ALL (k, i) zones (k=1..{n-3}, i=1..{n})...")

    results = []
    n_killed_by_higher_k = 0
    survives_all = []
    minimal_kills = {}
    t0 = time.time()
    for idx, ms in enumerate(step1_survivors_nontri):
        # Find any zone that kills this multiset
        zones = killing_zones(ms, n)
        if zones:
            # M was step-1 survivor, but some zone kills it. Note the minimal-k zone.
            min_k = min(k for (k, i) in zones)
            min_ki = min((k, i) for (k, i) in zones)
            minimal_kills[idx] = min_ki
            if any(k >= 2 for (k, _) in zones):
                n_killed_by_higher_k += 1
            results.append({'M': ms, 'killed_by': zones, 'min_killing_ki': min_ki})
        else:
            survives_all.append(ms)
            results.append({'M': ms, 'killed_by': [], 'min_killing_ki': None})
    dt = time.time() - t0
    print(f"  done in {dt:.1f}s")
    print(f"  step-1 survivors killed by some higher-k zone: "
          f"{n_killed_by_higher_k} / {len(step1_survivors_nontri)}")
    print(f"  step-1 survivors that ALSO survive every other zone: "
          f"{len(survives_all)} (expected 0 for World A)")
    return {
        'survives_all': survives_all,
        'n_killed_by_higher_k': n_killed_by_higher_k,
        'total_tested': len(step1_survivors_nontri),
        'per_survivor': results,
        'elapsed_seconds': dt,
    }


def phase1_for_n(n, only_orbit_reps=False, _orbit_rep_warn=True):
    """Full enumeration at this n. Returns survivors-of-everything set.
    only_orbit_reps: future hook for n=9 cyclic orbit reps; not implemented
    in this script (we run n=9 full)."""
    if only_orbit_reps and _orbit_rep_warn:
        print(f"  (orbit-rep enumeration not implemented; running full)")
    print(f"\n[Phase 1 at n={n}]  enumerating ALL multisets of size {n-3} "
          f"({len(list(combinations_with_replacement(all_chords(n), n-3)))} total)...")
    all_chord_list = all_chords(n)
    triangs        = []
    survives_all   = []
    killed_by_some = []
    minimal_kills_for_killed = []
    t0 = time.time()
    count = 0
    for ms in combinations_with_replacement(all_chord_list, n - 3):
        count += 1
        if is_triangulation(ms, n):
            triangs.append(ms)
            continue
        killed, ki = killed_by_some_zone(ms, n)
        if killed:
            killed_by_some.append(ms)
            minimal_kills_for_killed.append(ki)
        else:
            survives_all.append(ms)
        if count % 5000 == 0:
            print(f"  ... {count} processed ({time.time()-t0:.1f}s)")
    dt = time.time() - t0
    print(f"  done in {dt:.1f}s.  total={count}  triangulations={len(triangs)}  "
          f"killed-non-tri={len(killed_by_some)}  surviving-non-tri={len(survives_all)}")
    return {
        'n': n,
        'total_multisets':       count,
        'triangulations':        len(triangs),
        'killed_non_triangs':    len(killed_by_some),
        'surviving_non_triangs': len(survives_all),
        'surviving_list':        survives_all,
        'elapsed_seconds':       dt,
    }


# =========================================================================
# 7. Main
# =========================================================================

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
os.makedirs(OUT_DIR, exist_ok=True)


def write_n_findings(n, phase0_result, phase1_result, selfcheck_pass):
    md = []
    md.append(f"# Findings at n = {n}\n")
    md.append(f"- self-check (k=1 only reproduces step-1 counts): "
              f"**{'PASS' if selfcheck_pass else 'FAIL'}**")
    md.append(f"- Phase 0 (existing step-1 non-tri survivors tested vs all (k,i)):")
    md.append(f"    - total tested: {phase0_result['total_tested']}")
    md.append(f"    - killed by some higher-k zone: {phase0_result['n_killed_by_higher_k']}")
    md.append(f"    - **survives every (k,i)**: {len(phase0_result['survives_all'])}")
    md.append(f"    - elapsed: {phase0_result['elapsed_seconds']:.1f} s")
    if phase0_result['survives_all']:
        md.append(f"\nMultisets in Phase 0 that survive every (k,i):")
        for ms in phase0_result['survives_all']:
            md.append(f"  - {ms}")
    md.append(f"\n- Phase 1 (full enumeration):")
    md.append(f"    - total multisets: {phase1_result['total_multisets']}")
    md.append(f"    - triangulations: {phase1_result['triangulations']}")
    md.append(f"    - non-tri killed by some (k,i): {phase1_result['killed_non_triangs']}")
    md.append(f"    - **non-tri surviving every (k,i)**: "
              f"{phase1_result['surviving_non_triangs']}")
    md.append(f"    - elapsed: {phase1_result['elapsed_seconds']:.1f} s")
    if phase1_result['surviving_list']:
        md.append(f"\nNon-triangulations that survive every (k,i):")
        for ms in phase1_result['surviving_list'][:30]:
            md.append(f"  - {ms}")
        if len(phase1_result['surviving_list']) > 30:
            md.append(f"  (... {len(phase1_result['surviving_list'])-30} more)")
    return "\n".join(md) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ns', type=int, nargs='+', default=[7, 8],
                    help='which n values to run (default 7 8; add 9 for full but slow)')
    ap.add_argument('--skip-selfcheck', action='store_true')
    ap.add_argument('--skip-phase1', action='store_true',
                    help='Run Phase 0 only (Phase 1 full enumeration is slow at n=9)')
    args = ap.parse_args()

    expected_step1 = {5: 0, 6: 0, 7: 7, 8: 100, 9: 1011}
    overall = {}

    for n in args.ns:
        print(f"\n{'='*70}\n  n = {n}\n{'='*70}")
        # Self-check at n=7 and n=8 (n=9 is slow at full but we still
        # cross-check via step-1 survivors load below)
        if not args.skip_selfcheck:
            t0 = time.time()
            step1_nontri_survivors = selfcheck_step1(n, expected_step1.get(n, 0))
            step1_nontri_survivors = [ms for ms in step1_nontri_survivors
                                      if not is_triangulation(ms, n)]
            print(f"  (self-check took {time.time()-t0:.1f}s)")
        else:
            print(f"  (skipping self-check; recovering step-1 non-tri survivors)")
            step1_nontri_survivors = []
            for ms in combinations_with_replacement(all_chords(n), n - 3):
                if is_triangulation(ms, n):
                    continue
                if not any(killable_at_kzone(ms, n, 1, i) for i in range(1, n + 1)):
                    step1_nontri_survivors.append(ms)

        # Phase 0: take step-1 non-tri survivors, hit all zones
        p0 = phase0_for_n(n, step1_nontri_survivors)

        # Phase 1: full enumeration (skip if requested)
        if args.skip_phase1:
            print(f"\n[Phase 1 at n={n}]  SKIPPED (--skip-phase1)")
            p1 = {'n': n, 'total_multisets': None, 'triangulations': None,
                  'killed_non_triangs': None, 'surviving_non_triangs': None,
                  'surviving_list': [], 'elapsed_seconds': 0.0}
        else:
            p1 = phase1_for_n(n)

        # Save
        # Save with per-survivor minimal-(k,i) kill data
        per_survivor_serial = []
        for r in p0['per_survivor']:
            ms_list = [list(c) for c in r['M']]
            killed_by = [[k, i] for (k, i) in r['killed_by']]
            mki = list(r['min_killing_ki']) if r['min_killing_ki'] is not None else None
            per_survivor_serial.append({
                'M': ms_list,
                'killed_by': killed_by,
                'min_killing_ki': mki,
                'is_triangulation': is_triangulation(r['M'], n),
            })
        with open(os.path.join(OUT_DIR, f'n{n}_phase0.json'), 'w') as f:
            json.dump({
                'survives_all': [list(map(list, ms)) for ms in p0['survives_all']],
                'n_killed_by_higher_k': p0['n_killed_by_higher_k'],
                'total_tested': p0['total_tested'],
                'elapsed_seconds': p0['elapsed_seconds'],
                'per_survivor': per_survivor_serial,
            }, f, indent=2)
        with open(os.path.join(OUT_DIR, f'n{n}_phase1.json'), 'w') as f:
            json.dump({
                'n': p1['n'],
                'total_multisets': p1['total_multisets'],
                'triangulations': p1['triangulations'],
                'killed_non_triangs': p1['killed_non_triangs'],
                'surviving_non_triangs': p1['surviving_non_triangs'],
                'surviving_list': [list(map(list, ms)) for ms in p1['surviving_list']],
                'elapsed_seconds': p1['elapsed_seconds'],
            }, f, indent=2)
        md = write_n_findings(n, p0, p1, selfcheck_pass=not args.skip_selfcheck)
        with open(os.path.join(OUT_DIR, f'n{n}_findings.md'), 'w') as f:
            f.write(md)
        print(f"\n[n={n}] wrote n{n}_phase0.json, n{n}_phase1.json, n{n}_findings.md")
        overall[n] = {
            'phase0_survivors': len(p0['survives_all']),
            'phase1_surviving_non_triangs': p1['surviving_non_triangs'],
        }

    # Summary
    print(f"\n{'='*70}\n  VERDICT\n{'='*70}")
    world_a = all(v['phase1_surviving_non_triangs'] == 0 for v in overall.values())
    print("World A (no non-tri survives all leading kills): "
          f"{'YES' if world_a else 'NO'}")
    for n, v in overall.items():
        marker = "OK" if v['phase1_surviving_non_triangs'] == 0 else "PROBLEM"
        print(f"  n={n}: surviving non-triangulations = "
              f"{v['phase1_surviving_non_triangs']}  [{marker}]")
    summary = {'world_a': world_a, 'per_n': overall}
    with open(os.path.join(OUT_DIR, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nwrote outputs/summary.json")


if __name__ == "__main__":
    main()
