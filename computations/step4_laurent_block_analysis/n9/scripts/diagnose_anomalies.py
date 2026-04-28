"""
================================================================================
diagnose_anomalies.py  --  Deep structural diagnostic on n=9 step-1 survivors
                            that DON'T admit a depth-1 Laurent cascade
================================================================================

PURPOSE (BIG PICTURE)
---------------------
The depth-1 Laurent cascade verified the kill of every n=7 fish (7/7) and
every n=8 step-1 survivor (100/100, all 13 cyclic orbits).  At n=9 the
same depth-1 recipe failed to close orbit 22 after a 5705-second
exhaustive search, and four other orbits (25, 28, 34, 35) hit the
600 s per-rep timeout in the resumed run.  This script does the deep
diagnostic the user asked for:

  1.  Print the structural detail of M_rep at every zone Z:
      special, bare, all (companion, substitute) pairs at Z, which
      members of M_rep fall into each category, ell_Z(M_rep), and
      whether step-1 (K1)+(K2) hold (they should not, since M is a
      step-1 survivor).

  2.  Confirm the depth-1 search is exhaustive.  For every zone Z and
      every non-bare substitute U in M_rep at Z, compute the depth-1
      fingerprint  Y_U / (rest of M's free factors)  at order
      leading + 1, evaluate the resulting fingerprint equation, list
      cousins, and check step-1 killability of each cousin.

  3.  Try DEPTH-2 cascade at every zone.  For each zone with at least
      one non-bare substitute, try fingerprints
          Y_U^2 / (rest)               -- one substitute used twice in tail
          Y_U Y_{U'} / (rest)          -- two substitutes mixed in tail
      at order leading + 2.

  4.  If depth-2 fails, try DEPTH-3 and DEPTH-4 with the analogous
      fingerprints.

  5.  For each successful recipe, report (Z, depth, fingerprint,
      M_rep scalar, list of cousins and their step-1 kill zones).

  6.  Compare to the n=8 cascade structure: does M_rep have the same
      "missing-vertex forced ear" pattern as the n=7 fish, or does it
      need new structural ideas?

USAGE
-----
  python3 diagnose_anomalies.py 22                 # diagnose orbit 22
  python3 diagnose_anomalies.py 22 25 28 34 35     # diagnose multiple

Outputs go to `diagnose_orbit_<id>.md` per orbit and a combined
`anomaly_summary.md`.
"""

import json
import os
import sys
import time
from itertools import combinations_with_replacement, product

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.normpath(os.path.join(HERE, "..", "outputs"))
DIAG_DIR = os.path.normpath(os.path.join(HERE, "..", "diagnostics"))
N8_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "n8", "scripts"))
sys.path.insert(0, N8_DIR)

from cascade_kill_n8 import (  # noqa: E402
    all_chords, zone_structure, killable_at_zone_layer0, layer0_kill_zone,
    chord_symbol, laurent_coefficient, coefficient_of_monomial,
    candidate_multisets_for_fingerprint, format_multiset,
)
from analyze_recipes import vertices_of, missed_vertices  # noqa: E402


# ============================================================================
# 1.  STRUCTURAL CLASSIFICATION OF M AT EACH ZONE
# ============================================================================

def classify_chords_at_zone(M, zone_r, n):
    """
    Bucket the chords of M into {special, bare, substitutes (with companion
    info), free}.  Returns a dict per chord and a (n_special, n_bare,
    n_subs, n_free) tuple.

    LOGIC:  iterate M; for each chord, look up its role in
            zone_structure(zone_r, n).

    PHYSICS:  ell_Z(M) = n_special + n_bare + n_subs is the leading order
              in 1/X_S of M's expansion on Z.
    """
    special, bare, pairs = zone_structure(zone_r, n)
    sub_to_comp = {sub: comp for comp, sub in pairs}
    bucket = []
    n_special = n_bare = n_subs = n_free = 0
    for c in M:
        if c == special:
            bucket.append((c, "special", None))
            n_special += 1
        elif c == bare:
            bucket.append((c, "bare", None))
            n_bare += 1
        elif c in sub_to_comp:
            comp = sub_to_comp[c]
            bucket.append((c, "substitute", comp))
            n_subs += 1
        else:
            bucket.append((c, "free", None))
            n_free += 1
    return {
        "bucket": bucket,
        "n_special": n_special,
        "n_bare": n_bare,
        "n_subs": n_subs,
        "n_free": n_free,
        "ell_Z": n_special + n_bare + n_subs,
        "special": special,
        "bare": bare,
        "pairs": pairs,
    }


def step1_killable_diagnosis(M, n):
    """
    For each zone Z, return why M fails (K1)+(K2).

    LOGIC:  call classify_chords_at_zone; check (K1) (special or bare in M)
            and (K2) (some pair has both companion and substitute in M).

    PHYSICS:  certifies that M is genuinely a step-1 survivor.
    """
    out = []
    for r in range(1, n + 1):
        cls = classify_chords_at_zone(M, r, n)
        special = cls["special"]
        bare = cls["bare"]
        pairs = cls["pairs"]
        M_set = set(M)
        K1_fail_special = special in M_set
        K1_fail_bare = bare in M_set
        K2_fails = []
        for comp, sub in pairs:
            if comp in M_set and sub in M_set:
                K2_fails.append((comp, sub))
        killable = (not K1_fail_special and not K1_fail_bare and not K2_fails)
        out.append({
            "r": r,
            "special": special,
            "bare": bare,
            "ell_Z": cls["ell_Z"],
            "K1_fail_special": K1_fail_special,
            "K1_fail_bare": K1_fail_bare,
            "K2_fails": K2_fails,
            "killable": killable,
        })
    return out


# ============================================================================
# 2.  GENERIC DEPTH-d CASCADE SEARCH
# ============================================================================

def free_vars_at_zone(zone_r, n):
    """
    Return list of sympy symbols that are 'free' on zone Z_{zone_r, *}:
    every chord except special, bare, and the non-bare substitutes.
    """
    special, bare, pairs = zone_structure(zone_r, n)
    sub_set = {sub for _, sub in pairs}
    out = []
    for c in all_chords(n):
        if c == special or c == bare or c in sub_set:
            continue
        out.append(chord_symbol(c))
    return out


def build_depth_d_fingerprints(M, zone_r, n, depth):
    """
    Generate candidate fingerprints for a depth-`depth` cascade at zone Z.

    LOGIC:
      For each multiset of size `depth` chosen from the SUBSTITUTES of M
      at this zone (with replacement), build:
          numerator = product of companion symbols Y_{U_i} for chosen U_i's
          denominator = product of free symbols of M after removing one
              copy of each chosen substitute (with multiplicity).
      Yield (fp, chosen_substitutes) per candidate.

    PHYSICS:  at depth d, the Laurent expansion of 1/(Y_U - X_S) contributes
              a Y_U^d / X_S^{d+1} term, plus cross terms between different
              substitutes.  The fingerprint at order leading+d is therefore
              built from the chosen multiset of substitutes whose tails are
              "advanced" by one extra power each.  So the search space at
              depth d is "multisets of size d from the substitutes in M".
    """
    cls = classify_chords_at_zone(M, zone_r, n)
    sub_to_comp = {sub: comp for c, kind, comp in cls["bucket"]
                   if kind == "substitute"
                   for sub in [c]}
    # Substitutes IN M (with multiplicity, but we use the distinct list as
    # the "library" -- multiplicity just means same companion can be used).
    # Actually for the fingerprint we need to draw `depth` substitute
    # *occurrences* from M, but at the Laurent level we just care about
    # which substitute species; the multiplicity in M only constrains how
    # many independent occurrences are available.
    sub_chords = [c for c, kind, _ in cls["bucket"] if kind == "substitute"]
    if not sub_chords:
        return []

    # Multiset of size `depth` from sub_chords, but each draw cannot exceed
    # the multiplicity of that substitute in M.
    from collections import Counter
    sub_counter = Counter(sub_chords)
    distinct_subs = sorted(sub_counter.keys())

    def gen_choices(remaining_depth, idx):
        """Yield list of (chord, count) such that sum of counts == depth
        and each count <= sub_counter[chord]."""
        if remaining_depth == 0:
            yield []
            return
        if idx == len(distinct_subs):
            return
        max_take = min(sub_counter[distinct_subs[idx]], remaining_depth)
        for take in range(max_take + 1):
            for tail in gen_choices(remaining_depth - take, idx + 1):
                if take == 0:
                    yield tail
                else:
                    yield [(distinct_subs[idx], take)] + tail

    seen = set()
    out = []
    for choice in gen_choices(depth, 0):
        # Build the fingerprint.
        # Denominator: M's distinct chord factors (with multiplicities)
        # MINUS the chosen substitutes (with chosen counts), MINUS bares
        # and specials and unchosen substitutes (those don't show up in
        # the order-leading+depth fingerprint as free-variable factors).
        from collections import Counter as C2
        M_counter = C2(M)
        for sub_chord, take in choice:
            M_counter[sub_chord] -= take
            # keep only positive counts
        # Drop chord factors not contributing to free-variable denominator
        # at this Laurent order: special, bare, remaining substitutes.
        bare = cls["bare"]
        special = cls["special"]
        sub_set_local = {sub for _, sub in cls["pairs"]}
        denom_factors = []
        for c, cnt in M_counter.items():
            if cnt <= 0:
                continue
            if c == special or c == bare or c in sub_set_local:
                continue
            for _ in range(cnt):
                denom_factors.append(c)
        # Numerator: product of companion symbols Y_{U} raised to the
        # chosen count.
        numer_factors = []
        for sub_chord, take in choice:
            comp = None
            for cc, ss in cls["pairs"]:
                if ss == sub_chord:
                    comp = cc
                    break
            if comp is None:
                continue
            for _ in range(take):
                numer_factors.append(comp)
        # Build sympy fingerprint.
        fp = sp.Integer(1)
        for c in numer_factors:
            fp *= chord_symbol(c)
        for c in denom_factors:
            fp /= chord_symbol(c)
        # Canonical form for dedup.
        fp_simplified = sp.simplify(fp)
        key = str(fp_simplified)
        if key in seen:
            continue
        seen.add(key)
        out.append({
            "fp": fp_simplified,
            "fp_string": key,
            "chosen_substitutes": [(list(s), c) for s, c in choice],
            "numerator_chords": [list(c) for c in numer_factors],
            "denominator_chords": [list(c) for c in denom_factors],
        })
    return out


def fingerprint_equation_for_M(M, zone_r, n, target_order, fp,
                                pool_size_cap=None):
    """
    Build the fingerprint equation at (Z, target_order, fp).

    LOGIC:  candidate_multisets_for_fingerprint -> Laurent expansion ->
            scalar coefficient of fp.  Returns list of (multiset, scalar).

    PHYSICS:  exactly the depth-1 fingerprint_equation but generic in
              order; just iterates candidates and reads off coefficients.
    """
    pool = candidate_multisets_for_fingerprint(fp, n, n - 3)
    if pool_size_cap is not None and len(pool) > pool_size_cap:
        return None  # caller will detect and skip
    free_vars = free_vars_at_zone(zone_r, n)
    out = []
    for ms in pool:
        coeff_dict = laurent_coefficient(ms, zone_r, n, target_order)
        order_k_expr = coeff_dict.get(target_order, 0)
        if order_k_expr == 0:
            continue
        scalar = coefficient_of_monomial(order_k_expr, fp, free_vars)
        if scalar is None or scalar == 0:
            continue
        out.append((ms, scalar))
    return out


def try_depth_d_cascade(M, n, depth, pool_size_cap=10000, verbose=False):
    """
    For each zone Z, build all depth-`depth` fingerprints of M at Z, evaluate
    each, and check if M is in the equation with all cousins step-1 killable.

    Returns a list of recipe dicts (potentially empty).
    """
    out = []
    M_sorted = tuple(sorted(M))
    for zone_r in range(1, n + 1):
        cls = classify_chords_at_zone(M, zone_r, n)
        if cls["n_subs"] == 0:
            # No substitutes at this zone -> Laurent tails produce no
            # companion variables in the fingerprint.
            continue
        target_order = cls["ell_Z"] + depth
        fps = build_depth_d_fingerprints(M, zone_r, n, depth)
        if verbose:
            print(f"  zone {zone_r}: {len(fps)} fp candidates at order "
                  f"{target_order}")
        for fp_rec in fps:
            fp = fp_rec["fp"]
            t0 = time.time()
            eq = fingerprint_equation_for_M(M, zone_r, n, target_order, fp,
                                             pool_size_cap=pool_size_cap)
            if eq is None:
                if verbose:
                    print(f"    [skip] pool too large for "
                          f"{fp_rec['fp_string']}")
                continue
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
                cousin_kills.append({"multiset": [list(c) for c in m],
                                      "scalar": int(s),
                                      "step1_kill_zone": z0})
            elapsed = time.time() - t0
            recipe = {
                "zone_r": zone_r,
                "depth": depth,
                "order": target_order,
                "fingerprint": str(sp.simplify(fp)),
                "fp_chosen_substitutes": fp_rec["chosen_substitutes"],
                "fish_scalar": int(fish_scalar),
                "cousin_kills": cousin_kills,
                "all_cousins_step1_killable": all_dead,
                "elapsed_s": elapsed,
            }
            out.append(recipe)
            if verbose and all_dead:
                print(f"    [HIT] depth {depth} at zone {zone_r}: "
                      f"{fp_rec['fp_string']}")
    return out


# ============================================================================
# 3.  DRIVER
# ============================================================================

def diagnose_one(orbit_id, M_rep, n=9, max_depth=4, pool_cap=10000,
                 out_path=None):
    """
    Run the full diagnostic on one orbit representative.
    """
    lines = []
    def L(s):
        lines.append(s)
        print(s, flush=True)

    L(f"# Diagnostic for orbit {orbit_id}\n")
    L(f"M_rep = {format_multiset(tuple(tuple(c) for c in M_rep))}\n")
    M = tuple(tuple(c) for c in M_rep)
    V = sorted(vertices_of(M))
    Vm = sorted(missed_vertices(M, n))
    L(f"V_used = {V}\nV_missed = {Vm}\n")

    # Step-1 status across all zones.
    L("## Step-1 (K1)+(K2) status across all zones\n")
    diag = step1_killable_diagnosis(M, n)
    L("| zone | special | bare | ell_Z | K1 fail (special in M) | "
      "K1 fail (bare in M) | K2 fails (pairs) | killable? |")
    L("|---|---|---|---|---|---|---|---|")
    for d in diag:
        K2_str = ",".join(f"({c},{s})" for c, s in d["K2_fails"]) or "—"
        L(f"| Z_{{{d['r']},{(d['r']+2-1)%n+1}}} | {d['special']} | "
          f"{d['bare']} | {d['ell_Z']} | "
          f"{'YES' if d['K1_fail_special'] else 'no'} | "
          f"{'YES' if d['K1_fail_bare'] else 'no'} | "
          f"{K2_str} | {'killable' if d['killable'] else 'NOT killable'} |")
    L("")
    n_killable = sum(1 for d in diag if d["killable"])
    if n_killable == 0:
        L(f"Confirmed: M_rep is a genuine step-1 survivor (not killable "
          f"at any of the {n} zones).\n")
    else:
        L(f"WARNING: M_rep is killable at {n_killable} zones -- it should "
          f"NOT be a step-1 survivor!\n")

    # Try depth-1 first (sanity check), then 2, 3, 4.
    for depth in range(1, max_depth + 1):
        L(f"## Depth-{depth} cascade attempts\n")
        t0 = time.time()
        recipes = try_depth_d_cascade(M, n, depth, pool_size_cap=pool_cap,
                                       verbose=False)
        L(f"Tried depth-{depth} cascades across all zones in "
          f"{time.time()-t0:.1f} s; "
          f"found {len(recipes)} fingerprint equations containing M_rep "
          f"with nonzero scalar.")
        kill_recipes = [r for r in recipes if r["all_cousins_step1_killable"]]
        L(f"Of those, {len(kill_recipes)} have ALL cousins step-1 "
          f"killable (= valid kill recipe).\n")
        for r in kill_recipes:
            L(f"### Valid recipe @ depth {depth}, zone Z_{{{r['zone_r']},"
              f"{(r['zone_r']+2-1)%n+1}}}, order {r['order']}\n")
            L(f"- Fingerprint: `{r['fingerprint']}`")
            L(f"- Chosen substitutes (with multiplicity): "
              f"{r['fp_chosen_substitutes']}")
            L(f"- M_rep scalar: {r['fish_scalar']}")
            L(f"- Cousins ({len(r['cousin_kills'])}): all step-1 killable.")
            for c in r["cousin_kills"][:8]:
                L(f"  - {c['multiset']}: scalar={c['scalar']}, "
                  f"kill at Z_{{{c['step1_kill_zone']},*}}")
            if len(r["cousin_kills"]) > 8:
                L(f"  - ... ({len(r['cousin_kills'])-8} more)")
            L("")
        if kill_recipes:
            L(f"=== STOP: depth-{depth} closes orbit {orbit_id}. ===\n")
            if out_path:
                with open(out_path, "w") as f:
                    f.write("\n".join(lines) + "\n")
            return {
                "orbit_id": orbit_id,
                "depth_resolved": depth,
                "n_recipes": len(kill_recipes),
                "first_recipe": kill_recipes[0],
            }
        # Show a few non-closing recipes (cousins not all step-1 killable)
        # to see how close we are.
        partial = [r for r in recipes
                   if not r["all_cousins_step1_killable"]][:5]
        if partial:
            L(f"Top {len(partial)} partial fingerprint equations at depth "
              f"{depth} (cousins not all step-1 killable):\n")
            for r in partial:
                bad_cousins = [c for c in r["cousin_kills"]
                               if c["step1_kill_zone"] is None]
                L(f"- Z_{{{r['zone_r']},{(r['zone_r']+2-1)%n+1}}}, "
                  f"fp `{r['fingerprint']}`: M scalar={r['fish_scalar']}, "
                  f"{len(bad_cousins)} survivor-cousins out of "
                  f"{len(r['cousin_kills'])}")
                for c in bad_cousins[:3]:
                    L(f"  - survivor cousin: {c['multiset']} "
                      f"(scalar {c['scalar']})")
                L("")

    L(f"\n=== UNRESOLVED at depth <= {max_depth}. ===\n")
    if out_path:
        with open(out_path, "w") as f:
            f.write("\n".join(lines) + "\n")
    return {
        "orbit_id": orbit_id,
        "depth_resolved": None,
        "n_recipes": 0,
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: diagnose_anomalies.py ORBIT_ID [ORBIT_ID ...]")
        sys.exit(1)
    orbit_ids = [int(x) for x in sys.argv[1:]]

    with open(os.path.join(OUT_DIR, "orbits_n9.json")) as f:
        manifest = json.load(f)
    by_id = {o["orbit_id"]: o for o in manifest["orbits"]}

    results = []
    for oid in orbit_ids:
        if oid not in by_id:
            print(f"Orbit {oid} not in manifest.")
            continue
        out_path = os.path.join(DIAG_DIR, f"diagnose_orbit_{oid}.md")
        rec = diagnose_one(oid, by_id[oid]["representative"],
                            n=9, max_depth=4, pool_cap=10000,
                            out_path=out_path)
        results.append(rec)
        print(f"-> wrote {out_path}\n")

    # Combined summary.
    summary_lines = ["# Anomaly diagnostics summary\n"]
    summary_lines.append("| Orbit | Depth resolved? | # valid recipes |")
    summary_lines.append("|---:|---:|---:|")
    for r in results:
        d = r["depth_resolved"]
        summary_lines.append(f"| {r['orbit_id']} | "
                             f"{d if d is not None else '—'} | "
                             f"{r['n_recipes']} |")
    with open(os.path.join(HERE, "..", "anomaly_summary.md"), "w") as f:
        f.write("\n".join(summary_lines) + "\n")
    print("Wrote anomaly_summary.md.")


if __name__ == "__main__":
    main()
