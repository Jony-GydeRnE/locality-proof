"""
================================================================================
find_bridge_equations.py  --  Search for triangulation-only bridge equations
                                across Step-2 components at n=9
================================================================================

PURPOSE (BIG PICTURE)
---------------------
Triangulations are LOCAL terms; they are *not* meant to be killed. The
unitarity claim is that all 49 triangulation coefficients $c_T$ collapse
to a single common scalar $c$. Step-2 alone equates triangulations
within each of 16 components but not across; the cluster matrix at
depth-1 bridges 12 of the 16 components to each other; we are looking
for additional 1-zero equations that bridge the remaining 4 (the
"fan-class" components 1, 3, 7, 14).

A *bridge equation* of the form we want has the shape, after every
already-killable non-local coefficient is set to zero,
       sum_j  beta_j * c_{T_j}  =  0
where the $T_j$'s span at least TWO distinct Step-2 components, AT
LEAST ONE of which is in the untouched set {1, 3, 7, 14}. This forces
an equality of $c$-values across the two components when combined
with Step-2's within-component equating.

Every n=9 step-1 non-tri survivor is killed by some mechanism (Step-1
directly, single-orbit cascade, or block-rule cluster kill), so EVERY
non-local cousin reduces to 0. Hence every depth-1 or depth-2
fingerprint equation reduces to a pure triangulation relation. The
question is whether that relation is *trivial* (within one Step-2
component, hence already known via bare/special swap) or *bridging*
(spans multiple components).

USAGE
-----
  python3 find_bridge_equations.py [--depth 2] [--time-budget 3600]

OUTPUTS
-------
  outputs/n9_bridge_equations_depth<d>.json
      All bridge equations found at depth d.
  outputs/n9_bridges_summary.txt
      Final tally + which untouched components got bridged.
"""

import json
import os
import sys
import time
from collections import defaultdict, Counter

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
N9_OUT = os.path.normpath(os.path.join(
    HERE, "..", "..", "step4_laurent_block_analysis", "n9", "outputs"))
N8_DIR = os.path.normpath(os.path.join(
    HERE, "..", "..", "step4_laurent_block_analysis", "n8", "scripts"))
N9_SCRIPTS = os.path.normpath(os.path.join(
    HERE, "..", "..", "step4_laurent_block_analysis", "n9", "scripts"))
sys.path.insert(0, N8_DIR)
sys.path.insert(0, N9_SCRIPTS)

from cascade_kill_n8 import (  # noqa: E402
    layer0_kill_zone, is_triangulation, format_multiset, laurent_coefficient,
)
from analyze_recipes import canonical_orbit_rep  # noqa: E402
from cluster_analysis import depth1_fingerprint_equations  # noqa: E402
from diagnose_anomalies import (  # noqa: E402
    classify_chords_at_zone, build_depth_d_fingerprints,
    fingerprint_equation_for_M,
)

# Memoize Laurent expansions (single-process speedup).
import cascade_kill_n8 as _ckn8  # noqa: E402
import diagnose_anomalies as _diag  # noqa: E402
_LAURENT_CACHE = {}
_orig_laurent = _ckn8.laurent_coefficient


def _cached_laurent(monomial_chords, zone_r, n, max_order):
    key = (tuple(sorted(monomial_chords)), zone_r, n, max_order)
    if key not in _LAURENT_CACHE:
        _LAURENT_CACHE[key] = _orig_laurent(monomial_chords, zone_r, n,
                                             max_order)
    return _LAURENT_CACHE[key]


_ckn8.laurent_coefficient = _cached_laurent
_diag.laurent_coefficient = _cached_laurent

OUT_DIR = os.path.join(HERE, "outputs")


# ============================================================================
# Helpers
# ============================================================================

def load_orbit_manifest():
    """Map canonical-rep -> step-1 survivor orbit_id, and id -> rec."""
    with open(os.path.join(N9_OUT, "orbits_n9.json")) as f:
        d = json.load(f)
    rep_to_id = {}
    id_to_rec = {}
    for o in d["orbits"]:
        rep = tuple(tuple(c) for c in o["representative"])
        rep_to_id[rep] = o["orbit_id"]
        id_to_rec[o["orbit_id"]] = {"rep": rep, "size": o["size"]}
    return rep_to_id, id_to_rec


def load_step2_components():
    """Build orbit_id -> Step-2 component_id mapping for triangulations."""
    with open(os.path.join(HERE, "outputs",
                           "n9_triangulation_orbits.json")) as f:
        tri_data = json.load(f)
    rep_to_tri_id = {}
    for o in tri_data["orbits"]:
        rep = tuple(tuple(c) for c in o["representative"])
        rep_to_tri_id[rep] = o["orbit_id"]

    with open(os.path.join(HERE, "outputs",
                           "n9_step2_bare_swap_pairs.json")) as f:
        swap = json.load(f)
    edges = [(e["orbit_pair"][0], e["orbit_pair"][1])
             for e in swap["orbit_edges"]]
    parent = {oid: oid for oid in rep_to_tri_id.values()}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    # Assign sequential component IDs by min orbit ID per component.
    by_root = defaultdict(list)
    for oid in rep_to_tri_id.values():
        by_root[find(oid)].append(oid)
    sorted_roots = sorted(by_root.keys(), key=lambda r: min(by_root[r]))
    cid_map = {r: i + 1 for i, r in enumerate(sorted_roots)}
    rep_to_component = {
        rep: cid_map[find(rep_to_tri_id[rep])]
        for rep in rep_to_tri_id
    }
    return rep_to_component


# ============================================================================
# CORE: classify cousin in a fingerprint equation
# ============================================================================

def classify_cousin(ms, n, rep_to_step1_id, rep_to_tri_component):
    """
    Classify a cousin multiset.

    LOGIC:  if step-1 killable -> "step1_killed".
            elif is_triangulation -> "triangulation" (with component id).
            elif canonical rep is a step-1 non-tri survivor orbit
                -> "non_tri_survivor" (with orbit id; will be killed by
                its own cascade or by the cluster block rule).
            else -> "unknown" (shouldn't happen if classification is
                exhaustive).

    PHYSICS:  every n=9 multiset contributing to a fingerprint equation
              is one of: step-1 killable (a_M = 0 by leading-Laurent),
              a step-1 survivor non-triangulation (= 0 by cascade), a
              triangulation (NOT killed -- contributes c_T to the
              reduced equation), or "unknown" (= a chord multiset whose
              canonical rep isn't in our step-1-survivor manifest, which
              would indicate a manifest gap).
    """
    if layer0_kill_zone(list(ms), n) is not None:
        return ("step1_killed", None, None)
    if is_triangulation(list(ms), n):
        rep = canonical_orbit_rep(ms, n)
        cid = rep_to_tri_component.get(rep)
        if cid is None:
            return ("triangulation_unknown_component", rep, None)
        return ("triangulation", rep, cid)
    rep = canonical_orbit_rep(ms, n)
    oid = rep_to_step1_id.get(rep)
    if oid is not None:
        return ("non_tri_survivor", rep, oid)
    return ("unknown", rep, None)


def equation_to_reduced_triangulation_relation(
        eq, n, rep_to_step1_id, rep_to_tri_component):
    """
    Take a fingerprint equation (list of (multiset, scalar)) and reduce
    it by killing all non-triangulation cousins.

    LOGIC:  for each (ms, scalar) in eq, classify; if non-triangulation
            (step1_killed or non_tri_survivor), drop -- it contributes
            0; if triangulation, accumulate scalar by Step-2 component
            id (since within-component triangulations are equated).

    PHYSICS:  produces  sum_k  (sum-of-scalars-in-component-k) * c_k = 0
              where c_k is the within-component scalar (well-defined up
              to sign). If at least 2 distinct components have nonzero
              accumulated scalar, that's a BRIDGE.

    Returns:  dict component_id -> accumulated_scalar  (only nonzero).
              Plus a flag if any cousin was "unknown" (manifest gap).
    """
    comp_scalar = defaultdict(int)
    triangulation_terms = []
    has_unknown = False
    for ms, scalar in eq:
        kind, rep, info = classify_cousin(ms, n, rep_to_step1_id,
                                           rep_to_tri_component)
        if kind == "step1_killed":
            continue
        if kind == "non_tri_survivor":
            # Will be killed by cascade or block rule. Drop.
            continue
        if kind == "triangulation":
            cid = info
            comp_scalar[cid] += int(scalar)
            triangulation_terms.append({
                "rep": [list(c) for c in rep],
                "scalar": int(scalar),
                "step2_component": cid,
            })
            continue
        if kind == "triangulation_unknown_component":
            has_unknown = True
            continue
        if kind == "unknown":
            has_unknown = True
            continue
    # Drop zero-scalar components.
    comp_scalar = {k: v for k, v in comp_scalar.items() if v != 0}
    return {
        "component_scalars": dict(comp_scalar),
        "n_components": len(comp_scalar),
        "triangulation_terms": triangulation_terms,
        "has_unknown_cousin": has_unknown,
    }


# ============================================================================
# Driver
# ============================================================================

def find_bridges(depth, time_budget_s, target_untouched, log):
    """
    Enumerate fingerprint equations at given depth, reduce each by
    setting non-tri cousins to 0, and look for bridges that include
    at least one component in target_untouched.

    LOGIC:  iterate every step-1 survivor orbit (cluster + non-cluster);
            for each, build depth-d fingerprints across all zones;
            evaluate; reduce; if at least 2 components nonzero AND at
            least one in target_untouched, record bridge.

    PHYSICS:  each bridge equation forces an equality of c-values
              across two Step-2 components.  The target untouched
              components are exactly those we want to bridge.
    """
    rep_to_step1_id, id_to_rec = load_orbit_manifest()
    rep_to_tri_component = load_step2_components()

    log(f"At depth {depth}: scanning {len(id_to_rec)} step-1 non-tri "
        f"orbits...")
    bridges = []
    untouched_reached = set()
    components_seen = set()
    n = 9
    t0 = time.time()
    for i, oid in enumerate(sorted(id_to_rec.keys())):
        if time.time() - t0 > time_budget_s:
            log(f"  TIME BUDGET EXCEEDED at orbit {i+1}/{len(id_to_rec)}; "
                f"stopping.")
            break
        M = id_to_rec[oid]["rep"]
        n_eqs = 0
        n_bridges_orbit = 0
        for zone_r in range(1, n + 1):
            cls = classify_chords_at_zone(M, zone_r, n)
            if cls["n_subs"] == 0 and cls["n_bare"] == 0:
                # No bare or substitute -> no Laurent expansion at this
                # zone gives us anything new.
                continue
            if depth == 1:
                fps = []
                for fp_rec in depth1_fingerprint_equations(M, n):
                    if fp_rec["zone_r"] == zone_r:
                        fps.append({
                            "fp": fp_rec["fp_expr"],
                            "fp_string": fp_rec["fingerprint"],
                        })
            else:
                fps = build_depth_d_fingerprints(M, zone_r, n, depth)
            target_order = cls["ell_Z"] + depth
            for fp_rec in fps:
                if depth == 1:
                    fp = fp_rec["fp"]
                    fp_string = fp_rec["fp_string"]
                else:
                    fp = fp_rec["fp"]
                    fp_string = fp_rec["fp_string"]
                eq = fingerprint_equation_for_M(
                    M, zone_r, n, target_order, fp,
                    pool_size_cap=8000)
                if eq is None or not eq:
                    continue
                # Reduce.
                red = equation_to_reduced_triangulation_relation(
                    eq, n, rep_to_step1_id, rep_to_tri_component)
                if red["has_unknown_cousin"]:
                    # Manifest gap -- skip; we can't reduce safely.
                    continue
                comp_scalars = red["component_scalars"]
                if red["n_components"] < 2:
                    continue
                touched_in_target = set(comp_scalars.keys()) & target_untouched
                components_seen.update(comp_scalars.keys())
                if not touched_in_target:
                    continue
                # BRIDGE
                bridges.append({
                    "from_orbit": oid,
                    "from_orbit_rep": [list(c) for c in M],
                    "zone": zone_r,
                    "depth": depth,
                    "Laurent_order": target_order,
                    "fingerprint": fp_string,
                    "component_scalars": comp_scalars,
                    "triangulation_terms": red["triangulation_terms"],
                    "untouched_in_eq": sorted(touched_in_target),
                })
                untouched_reached.update(touched_in_target)
                n_bridges_orbit += 1
                n_eqs += 1
        elapsed = time.time() - t0
        if (i + 1) % 5 == 0 or n_bridges_orbit > 0:
            log(f"  [{i+1}/{len(id_to_rec)}] orbit {oid}: "
                f"{n_bridges_orbit} bridges (cumul "
                f"{len(bridges)}; untouched reached "
                f"{sorted(untouched_reached)}; elapsed "
                f"{elapsed:.0f} s)")
    log(f"  done.  total bridges = {len(bridges)};  "
        f"untouched reached = {sorted(untouched_reached)}")
    return bridges, untouched_reached, components_seen


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--depth", type=int, default=2)
    p.add_argument("--time-budget", type=int, default=3600)
    args = p.parse_args()

    log_lines = []

    def log(s):
        log_lines.append(s)
        print(s, flush=True)

    log("=" * 78)
    log(f"Searching for triangulation-only BRIDGE equations at "
        f"depth {args.depth}.")
    log(f"Target: bridge any of the 4 untouched components "
        f"{{1, 3, 7, 14}} to another component.")
    log("=" * 78)

    target_untouched = {1, 3, 7, 14}
    bridges, untouched_reached, components_seen = find_bridges(
        args.depth, args.time_budget, target_untouched, log)

    out_json = os.path.join(OUT_DIR,
                             f"n9_bridge_equations_depth{args.depth}.json")
    with open(out_json, "w") as f:
        json.dump({
            "depth": args.depth,
            "n_bridges_found": len(bridges),
            "untouched_components_reached": sorted(untouched_reached),
            "untouched_components_remaining":
                sorted(target_untouched - untouched_reached),
            "bridges": bridges,
        }, f, indent=2)
    log(f"\nWrote {out_json}")

    out_summary = os.path.join(OUT_DIR, "n9_bridges_summary.txt")
    with open(out_summary, "w") as f:
        f.write("\n".join(log_lines) + "\n")
        f.write(f"\nDepth {args.depth} results:\n")
        f.write(f"  bridges found: {len(bridges)}\n")
        f.write(f"  untouched reached: {sorted(untouched_reached)}\n")
        f.write(f"  untouched still remaining: "
                f"{sorted(target_untouched - untouched_reached)}\n")
    log(f"Summary written to {out_summary}")


if __name__ == "__main__":
    main()
