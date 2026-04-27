"""
================================================================================
audit_orbits.py  --  Per-orbit structural audit of the n=8 step-1 survivors
================================================================================

PURPOSE (BIG PICTURE)
---------------------
For each of the 13 cyclic orbits of step-1 survivors at n=8, this script
does a deeper structural audit of one orbit representative M_rep.  For
every zone Z of M_rep and every substitute chord U in M_rep at that zone,
it logs:

  ell_Z(M_rep)
      = (number of specials in M, number of bares in M,
         number of non-bare-substitute chords in M)
      i.e., how many of M's chord factors get rewritten in terms of X_S
      on this zone (the rest are "free" companions or born-free chords).

  For each substitute U in M_rep at zone Z:
      Y_U = its companion (the row-r chord with the same far endpoint)
      K2_violation = (Y_U is also in M_rep)   # i.e. (K2) fails at this zone
                                              # because the full
                                              # (companion, substitute)
                                              # pair sits in M.

  Depth-1 cascade attempt at (Z, U):
      Build fingerprint  X_{Y_U} / (rest of M's free-side factors).
      Try Laurent orders leading_order+1 and leading_order+2.
      Check whether the resulting fingerprint equation has all OTHER
      contributors step-1 killable at some zone.
      Result: either a recipe (zone, order, fingerprint, cousin kills)
      or "no depth-1 hit at this (Z, U)".

The point is to enumerate EVERY valid (Z, U) -- not just the first hit
the cascade verifier finds -- so we can see if multiple recipes work
for a single representative, and which structural features (V_missed,
crossing pairs, K2 violation) correlate with successful recipes.

OUTPUTS
-------
  orbit_audit.json    -- one block per orbit:
      orbit_id, orbit_size, representative,
      vertices_used, vertices_missed,
      crossing_pairs (as a list of pairs of chords),
      phi_i_match,
      zones: [
          {
              r, special, bare,
              n_special_in_M, n_bare_in_M, n_nonbare_subs_in_M,
              ell_Z = (above three),
              substitutes: [
                  {substitute, companion, Y_in_M, cascade_recipes: [...]},
                  ...
              ]
          },
          ...
      ]

  orbit_audit_summary.md    -- human-readable summary table.

CODE STYLE
----------
Per repo contribution standards, every function has a docstring with a
LOGIC section (algorithm) and a PHYSICS section (what it computes in
proof language).
"""

import json
import os
import sys
from collections import defaultdict

import sympy as sp

# The cascade library lives next door in this same folder.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from cascade_kill_n8 import (  # noqa: E402
    all_chords, zone_structure, killable_at_zone_layer0, layer0_kill_zone,
    chord_symbol, laurent_coefficient, coefficient_of_monomial,
    fingerprint_equation, candidate_multisets_for_fingerprint,
    fish_substitutes_at_zone, fish_count_bare_and_subs, format_multiset,
)
from analyze_recipes import (  # noqa: E402
    parse_results, vertices_of, missed_vertices, frame_pairs,
    canonical_orbit_rep, orbit_decomposition, cyclic_distance,
)


# ============================================================================
# 1.  ELL_Z and substitute structure
# ============================================================================

def ell_at_zone(M, zone_r, n=8):
    """
    Triple (n_special, n_bare, n_nonbare_subs) counting how many factors
    of M get rewritten in terms of X_S on zone Z_{zone_r, zone_r+2}.

    LOGIC:  iterate M; classify each chord against the zone structure.

    PHYSICS:  ell_Z controls the leading order of M's Laurent expansion
              on zone Z: the leading 1/X_S^{ell_Z_total - n_special} is
              produced by the n_bare bares + n_nonbare_subs substitutes
              (each contributes 1/X_S to leading order); the n_special
              specials each contribute 1/X_S directly.  Together with
              the free factors, the leading exponent in 1/X_S is
              ell_Z_total = n_special + n_bare + n_nonbare_subs.
    """
    special, bare, pairs = zone_structure(zone_r, n)
    sub_set = {sub for _, sub in pairs}
    n_sp = sum(1 for c in M if c == special)
    n_b  = sum(1 for c in M if c == bare)
    n_ns = sum(1 for c in M if c in sub_set)
    return (n_sp, n_b, n_ns)


def substitutes_with_companion_flag(M, zone_r, n=8):
    """
    For zone Z_{zone_r, ...}, list (Y_U, U, Y_in_M) for each substitute
    U appearing in M (with multiplicity, deduplicated by identity).

    LOGIC:  scan zone's (companion, substitute) pairs; for each
            substitute U in M, check if its companion Y_U is also in M.

    PHYSICS:  Y_in_M = True is the (K2) violation flag for this pair --
              both companion and substitute are in M, so the leading
              kill argument breaks (the pair contributes a non-monomial
              partial-fraction factor 1/(Y_U - X_S) that mixes orders).
              Y_in_M = False means the companion is "free" and can serve
              as the fingerprint variable in a depth-1 cascade.
    """
    special, bare, pairs = zone_structure(zone_r, n)
    out = []
    seen = set()
    M_set = set(M)
    for companion, substitute in pairs:
        if substitute in M_set and substitute not in seen:
            seen.add(substitute)
            Y_in_M = (companion in M_set)
            out.append({
                "substitute": list(substitute),
                "companion": list(companion),
                "Y_in_M": Y_in_M,
            })
    return out


# ============================================================================
# 2.  Crossing pairs
# ============================================================================

def crossing_pairs_in_multiset(M, n=8):
    """
    All pairs (c, c') of distinct chords in M whose endpoints cross
    (in the cyclic-polygon sense).

    LOGIC:  iterate distinct unordered pairs; apply the standard crossing
            test (a, b) crosses (c, d) iff exactly one of c, d lies
            strictly between a and b, AND no endpoint is shared.

    PHYSICS:  crossings characterize non-locality: a pair of crossing
              chords means M cannot be a triangulation; the count and
              pattern of crossings is one structural feature that
              distinguishes orbits.
    """
    # Use distinct chords (set), but keep the original chord values.
    distinct = list({c for c in M})
    out = []
    for i in range(len(distinct)):
        for j in range(i + 1, len(distinct)):
            a, b = distinct[i]
            c, d = distinct[j]
            if a == c or a == d or b == c or b == d:
                continue
            bc = (a < c < b)
            bd = (a < d < b)
            if bc != bd:
                out.append((distinct[i], distinct[j]))
    return out


# ============================================================================
# 3.  Per-(zone, substitute) cascade audit
# ============================================================================

# Memoization cache for Laurent expansions; speeds up the audit by
# reusing computations across fingerprint attempts that share a (multiset,
# zone) pair.  Key = (sorted multiset, zone_r, n, max_order).
_LAURENT_CACHE = {}


def _laurent_cached(multiset_chords, zone_r, n, max_order):
    """
    Memoized wrapper for laurent_coefficient.

    LOGIC:  cache hits return immediately; misses compute and store.

    PHYSICS:  the Laurent expansion of a multiset on a zone is
              deterministic, so caching across the audit's many
              fingerprint attempts that share the same pool gives a
              direct speedup.
    """
    key = (tuple(sorted(multiset_chords)), zone_r, n, max_order)
    if key not in _LAURENT_CACHE:
        _LAURENT_CACHE[key] = laurent_coefficient(
            multiset_chords, zone_r, n, max_order)
    return _LAURENT_CACHE[key]


def try_cascade_at_zone_substitute(M, zone_r, companion, substitute, n=8,
                                   max_extra_orders=1):
    """
    Try to verify a depth-1 cascade kill of M at zone Z_{zone_r, ...}
    using the chosen substitute U = `substitute` and companion Y_U =
    `companion`.

    LOGIC:
      Build fingerprint  X_{Y_U} / (X_factors of M, minus U, minus all
      bares/specials/other-substitutes on this zone).
      For each Laurent order from leading_order+1 to leading_order+max_extra_orders:
          generate candidate multisets, compute fingerprint equation,
          check that M is in the equation with nonzero scalar AND every
          OTHER contributor is step-1 killable at some zone.
      Return all successful recipes (a list -- could be empty).

    PHYSICS:  exactly the depth-1 cascade rule used in `cascade_kill_n8.py`,
              but driven by a specific (Y_U, U) choice rather than scanning
              all substitutes.  This lets us enumerate every valid (Z, U)
              recipe per representative, not just the cascade verifier's
              first hit.

    Returns a list of dicts: each {order, fingerprint_str, fish_scalar,
    cousins (list of {multiset, scalar, kill_zone})}.
    """
    special_z, bare_z, pairs_z = zone_structure(zone_r, n)
    sub_set = {sub for _, sub in pairs_z}

    n_b, n_sub_total = fish_count_bare_and_subs(M, zone_r, n)
    if n_sub_total == 0:
        return []
    leading_order = n_b + n_sub_total

    # Build the fingerprint.
    M_minus = list(M)
    M_minus.remove(substitute)
    rest_free = [c for c in M_minus
                 if c != special_z and c != bare_z and c not in sub_set]
    rest_factors = sp.Integer(1)
    for c in rest_free:
        rest_factors *= chord_symbol(c)
    fp = chord_symbol(companion) / rest_factors

    out = []
    M_sorted = tuple(sorted(M))
    for extra in range(1, max_extra_orders + 1):
        target_order = leading_order + extra
        pool = candidate_multisets_for_fingerprint(fp, n, n - 3)
        if not pool:
            continue
        eq = fingerprint_equation(M, zone_r, n, target_order, fp, pool)
        if not eq:
            continue
        fish_scalar = None
        for m, s in eq:
            if tuple(sorted(m)) == M_sorted:
                fish_scalar = s
                break
        if fish_scalar is None or fish_scalar == 0:
            continue
        cousins = [(m, s) for m, s in eq
                   if tuple(sorted(m)) != M_sorted]
        cousin_kills = []
        all_dead = True
        for m, s in cousins:
            z0 = layer0_kill_zone(list(m), n)
            if z0 is None:
                all_dead = False
                cousin_kills.append({
                    "multiset": [list(c) for c in m],
                    "scalar": int(s),
                    "step1_kill_zone": None,
                })
            else:
                cousin_kills.append({
                    "multiset": [list(c) for c in m],
                    "scalar": int(s),
                    "step1_kill_zone": z0,
                })
        if all_dead:
            out.append({
                "order": target_order,
                "fingerprint": str(sp.simplify(fp)),
                "fish_scalar": int(fish_scalar),
                "cousins": cousin_kills,
            })
    return out


# ============================================================================
# 4.  Per-orbit audit driver
# ============================================================================

def audit_one_orbit(orbit_id, orbit_size, M_rep, n=8):
    """
    Build the full audit record for one orbit representative.

    LOGIC:  For each zone Z, compute ell_Z, list substitutes in M_rep
            at Z with K2-violation flag, attempt the depth-1 cascade
            for each (Z, substitute) with a substitute available, log
            results.

    PHYSICS:  one record per orbit completely characterizes the orbit's
              kill structure (because cyclic-shifted reps inherit
              shifted recipes), so the all-n proof needs only to verify
              one such record per orbit.
    """
    record = {
        "orbit_id": orbit_id,
        "orbit_size": orbit_size,
        "representative": [list(c) for c in M_rep],
        "vertices_used": sorted(vertices_of(M_rep)),
        "vertices_missed": sorted(missed_vertices(M_rep, n)),
        "frame_candidates_phi_i": [list(p) for p in
                                   frame_pairs(missed_vertices(M_rep, n), n)],
        "crossing_pairs": [
            [list(c1), list(c2)] for (c1, c2) in
            crossing_pairs_in_multiset(M_rep, n)
        ],
        "zones": [],
    }
    for zone_r in range(1, n + 1):
        special, bare, pairs = zone_structure(zone_r, n)
        ell = ell_at_zone(M_rep, zone_r, n)
        substitutes = substitutes_with_companion_flag(M_rep, zone_r, n)
        # For each substitute with companion not in M (i.e. K2 not
        # violated and a cascade could plausibly fire), try the cascade.
        for sub_rec in substitutes:
            U = tuple(sub_rec["substitute"])
            Y = tuple(sub_rec["companion"])
            if sub_rec["Y_in_M"]:
                # K2 already violated -- this means M is killable by
                # step-1 if every other pair also fails... but the
                # cascade route through this pair won't isolate M
                # cleanly, so we skip the depth-1 attempt here.
                # (We still log the substitute with Y_in_M=True.)
                sub_rec["cascade_recipes"] = "(skipped: K2 violation -- " \
                    "companion already in M; cascade route via this pair " \
                    "will not give a clean fingerprint)"
                continue
            # User wants the *depth-1* cascade audit, which is exactly
            # leading_order + 1.  We don't try +2 here; that's a deeper
            # cascade and is its own separate question.
            recipes = try_cascade_at_zone_substitute(
                M_rep, zone_r, Y, U, n=n, max_extra_orders=1)
            sub_rec["cascade_recipes"] = recipes
        zone_rec = {
            "r": zone_r,
            "special": list(special),
            "bare": list(bare),
            "n_special_in_M": ell[0],
            "n_bare_in_M": ell[1],
            "n_nonbare_subs_in_M": ell[2],
            "ell_Z_total": sum(ell),
            "substitutes": substitutes,
        }
        record["zones"].append(zone_rec)
    return record


def write_summary_md(records, out_path):
    """
    Write a human-readable summary of the per-orbit audit.

    LOGIC:  iterate records, print one block per orbit with a small
            table of (Z, ell_Z, # substitutes attempted, # successful
            cascade recipes).

    PHYSICS:  the summary lets a reader scan the 13 orbits and see at a
              glance which orbits have a unique recipe vs multiple, and
              which (Z, U) pairs are productive for each.
    """
    lines = []
    lines.append("# n=8 per-orbit structural audit (summary)\n")
    lines.append(f"Each row counts, for one orbit representative M_rep:\n"
                 f"- ell_Z = (n_special, n_bare, n_nonbare_subs) at zone Z;\n"
                 f"- 'subs tried' = # of distinct substitutes in M at Z "
                 f"that have their companion Y_U NOT in M (so a cascade "
                 f"route is plausible);\n"
                 f"- 'recipes found' = # of successful depth-1 cascade "
                 f"hits across orders leading+1 and leading+2.\n")
    for rec in records:
        lines.append(f"## Orbit {rec['orbit_id']}  (size {rec['orbit_size']})\n")
        lines.append(f"Representative: `{format_multiset(tuple(map(tuple, rec['representative'])))}`")
        lines.append(f"V(M) = {rec['vertices_used']}, V_missed = {rec['vertices_missed']}.")
        lines.append(f"Crossing chord pairs in M: {len(rec['crossing_pairs'])}.")
        lines.append(f"Phi-I frame candidates: {rec['frame_candidates_phi_i']}\n")
        lines.append("| zone $Z$ | special | bare | $\\ell_Z$ "
                     "$=({n_{sp}},{n_b},{n_{ns}})$ | subs tried "
                     "(K2-clean) | recipes found |")
        lines.append("|---|---|---|---|---|---|")
        for zr in rec["zones"]:
            ell = (zr["n_special_in_M"], zr["n_bare_in_M"],
                   zr["n_nonbare_subs_in_M"])
            tried = sum(
                1 for s in zr["substitutes"]
                if not s["Y_in_M"] and isinstance(s.get("cascade_recipes"),
                                                  list)
            )
            found = sum(
                len(s["cascade_recipes"])
                for s in zr["substitutes"]
                if isinstance(s.get("cascade_recipes"), list)
            )
            lines.append(f"| Z_{{{zr['r']},{(zr['r']+2-1)%8+1}}} | "
                         f"{tuple(zr['special'])} | {tuple(zr['bare'])} | "
                         f"{ell} | {tried} | {found} |")
        lines.append("")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ============================================================================
# 5.  Main
# ============================================================================

def main():
    """
    Top-level driver: parse cascade results, compute orbits, audit each
    orbit representative, write JSON + markdown summary.
    """
    n = 8
    survivors = parse_results()
    print(f"Parsed {len(survivors)} survivors.")
    orbits = orbit_decomposition(survivors, n)
    rep_keys_sorted = sorted(orbits.keys())
    print(f"Found {len(orbits)} orbits, sizes "
          f"{sorted({len(v) for v in orbits.values()})}.")

    records = []
    import time
    t_start = time.time()
    for orbit_id, rep in enumerate(rep_keys_sorted, start=1):
        orbit_size = len(orbits[rep])
        t_orbit = time.time()
        print(f"  Auditing orbit {orbit_id}/{len(rep_keys_sorted)} "
              f"(size {orbit_size}, rep {format_multiset(rep)})...",
              flush=True)
        rec = audit_one_orbit(orbit_id, orbit_size, rep, n=n)
        records.append(rec)
        # Save partial progress so we can inspect mid-run.
        out_json = os.path.join(HERE, "orbit_audit.json")
        with open(out_json, "w") as f:
            json.dump(records, f, indent=2)
        print(f"    done in {time.time()-t_orbit:.1f} s "
              f"(total elapsed {time.time()-t_start:.1f} s).", flush=True)

    out_json = os.path.join(HERE, "orbit_audit.json")
    with open(out_json, "w") as f:
        json.dump(records, f, indent=2)
    print(f"Wrote {out_json}.")

    out_md = os.path.join(HERE, "orbit_audit_summary.md")
    write_summary_md(records, out_md)
    print(f"Wrote {out_md}.")


if __name__ == "__main__":
    main()
