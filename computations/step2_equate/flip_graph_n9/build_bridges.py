"""
================================================================================
build_bridges.py  --  Bridge the 4 untouched Step-2 components at n=9
================================================================================

PURPOSE (BIG PICTURE)
---------------------
The Step-2 flip graph on the 49 triangulation orbits at n=9 has 16
connected components.  The 90-orbit cluster matrix at depth-1 from
`../../step4_laurent_block_analysis/n9/` provides 32 external
triangulation columns; those 32 orbits land in 12 of the 16
Step-2 components, leaving 4 components untouched.

We want to bridge those 4 untouched components into the cluster's
reach using ONLY 1-zero hidden-zero constraints (Rodina's setup;
no k-zeros for k >= 2).

This script runs the bridging tests in order.  It STOPS once
16/16 components are touched.

  TEST 1.  Depth-1 fingerprint equations of the 23 non-cluster
           non-tri step-1 survivor orbits.  These orbits were
           killed by their own singleton cascades, but their
           OTHER (zone, substitute) fingerprint attempts may
           involve triangulation cousins from previously
           untouched components.

  TEST 2.  Depth-2 (Laurent order leading + 2) fingerprint
           equations of every cluster + non-cluster non-tri orbit.
           Whenever a triangulation cousin from an untouched
           component appears in a depth-2 equation, that's a bridge.

  TEST 3.  Pure-triangulation Laurent extractions.  At each cyclic
           zone, at low Laurent orders, find free-chord monomials
           whose only contributors are triangulations (no non-local
           cousins).  Each yields a triangulation-only relation.

  TEST 4.  Structural characterization of the (still-)untouched
           components: cyclic stabilizer, chord-length distribution,
           ear vertices, common features.

OUTPUTS
-------
  outputs/n9_bridge_test_summary.txt   -- per-test result + verdict
  outputs/n9_bridge_test1_equations.json -- new triangulation cousins from TEST 1
  outputs/n9_bridge_test2_equations.json -- (only written if TEST 1 fails)
  outputs/n9_bridge_test3_equations.json -- (only written if TEST 2 fails)
  outputs/n9_untouched_structural.md    -- (only written if all tests fail)

CODE STYLE
----------
Per repo standards: every function has LOGIC + PHYSICS docstring.
"""

import json
import os
import sys
import time
from collections import defaultdict

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
    layer0_kill_zone, is_triangulation, format_multiset,
    laurent_coefficient,
)
from analyze_recipes import canonical_orbit_rep  # noqa: E402
from cluster_analysis import depth1_fingerprint_equations  # noqa: E402
from diagnose_anomalies import (  # noqa: E402
    classify_chords_at_zone, build_depth_d_fingerprints,
    fingerprint_equation_for_M,
)

# Monkey-patch laurent_coefficient with a memoized version for speed.
import cascade_kill_n8 as _ckn8  # noqa: E402
import diagnose_anomalies as _diag  # noqa: E402
import cluster_analysis as _clu  # noqa: E402
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
_clu.laurent_coefficient = _cached_laurent

OUT_DIR = os.path.join(HERE, "outputs")


# ============================================================================
# Setup helpers
# ============================================================================

def load_orbit_manifest():
    """Map canonical-rep -> orbit_id and orbit_id -> rec for n=9 step-1
    non-tri survivor orbits."""
    with open(os.path.join(N9_OUT, "orbits_n9.json")) as f:
        d = json.load(f)
    rep_to_id = {}
    id_to_rec = {}
    for o in d["orbits"]:
        rep = tuple(tuple(c) for c in o["representative"])
        rep_to_id[rep] = o["orbit_id"]
        id_to_rec[o["orbit_id"]] = {"rep": rep, "size": o["size"]}
    return rep_to_id, id_to_rec


def load_cluster_orbits():
    """Return the set of 90 cluster orbit IDs."""
    with open(os.path.join(N9_OUT, "cluster_partial.json")) as f:
        d = json.load(f)
    return set(d["cluster_orbits"])


def load_step2_components():
    """
    Build orbit_id -> Step-2 component_id mapping for the 49
    triangulation orbits.

    LOGIC:  load orbits + edges from flip_graph_n9 outputs; run
            union-find; assign component IDs starting at 1.

    PHYSICS:  components are the equivalence classes of triangulation
              coefficients under Step-2 bare/special swaps.  An orbit's
              component ID is the invariant we want to track when
              bridging.
    """
    with open(os.path.join(HERE, "outputs",
                           "n9_triangulation_orbits.json")) as f:
        tri_data = json.load(f)
    rep_to_tri_id = {}
    tri_id_to_rep = {}
    for o in tri_data["orbits"]:
        rep = tuple(tuple(c) for c in o["representative"])
        rep_to_tri_id[rep] = o["orbit_id"]
        tri_id_to_rep[o["orbit_id"]] = rep

    with open(os.path.join(HERE, "outputs",
                           "n9_step2_bare_swap_pairs.json")) as f:
        swap = json.load(f)
    edges = [(e["orbit_pair"][0], e["orbit_pair"][1])
             for e in swap["orbit_edges"]]

    # union-find
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

    rep_to_component = {}
    component_root_to_id = {}
    next_cid = 1
    for rep, tri_id in rep_to_tri_id.items():
        root = find(tri_id)
        if root not in component_root_to_id:
            component_root_to_id[root] = next_cid
            next_cid += 1
        rep_to_component[rep] = component_root_to_id[root]

    return rep_to_component, rep_to_tri_id


def load_initial_touched_components():
    """
    Return the set of Step-2 component IDs already touched by the
    cluster matrix's 32 external triangulation orbits.
    """
    with open(os.path.join(N9_OUT, "cluster_partial.json")) as f:
        d = json.load(f)
    rep_to_component, _ = load_step2_components()
    touched = set()
    for ec in d["external_columns"]:
        rep = tuple(tuple(c) for c in ec["rep"])
        if is_triangulation(list(rep), 9):
            cid = rep_to_component.get(rep)
            if cid is not None:
                touched.add(cid)
    return touched


# ============================================================================
# TEST 1 -- non-cluster orbit depth-1 equations
# ============================================================================

def test1_non_cluster_bridges(rep_to_component, rep_to_id, id_to_rec,
                               cluster_ids, log):
    """
    Enumerate every depth-1 fingerprint equation for each of the 23
    non-cluster non-tri step-1 survivor orbits.  For each equation,
    list any triangulation cousin and which Step-2 component it lies
    in.  Aggregate the set of components freshly touched by these
    equations.

    LOGIC:  for each non-cluster orbit, call
            depth1_fingerprint_equations(M_rep, 9); filter cousins
            for is_triangulation; map to component IDs.

    PHYSICS:  the non-cluster orbits' depth-1 equations may couple to
              triangulation cousins from previously untouched Step-2
              components, providing the "missing bridge" we need.
    """
    n = 9
    non_cluster_ids = sorted(set(id_to_rec.keys()) - cluster_ids)
    log(f"TEST 1: enumerating depth-1 equations for "
        f"{len(non_cluster_ids)} non-cluster orbits...")

    new_touched = set()
    bridge_records = []
    t0 = time.time()
    for i, oid in enumerate(non_cluster_ids):
        M = id_to_rec[oid]["rep"]
        n_eqs = 0
        for eq in depth1_fingerprint_equations(M, n):
            tri_cousins = []
            for term in eq["equation"]:
                ms = tuple(tuple(c) for c in term["multiset"])
                if not is_triangulation(list(ms), n):
                    continue
                rep = canonical_orbit_rep(ms, n)
                cid = rep_to_component.get(rep)
                if cid is not None:
                    tri_cousins.append({
                        "tri_rep": [list(c) for c in rep],
                        "scalar": term["scalar"],
                        "step2_component": cid,
                    })
                    new_touched.add(cid)
            if tri_cousins:
                bridge_records.append({
                    "from_orbit": oid,
                    "from_orbit_rep": [list(c) for c in M],
                    "zone": eq["zone_r"],
                    "fingerprint": eq["fingerprint"],
                    "tri_cousins": tri_cousins,
                })
            n_eqs += 1
        log(f"  [{i+1}/{len(non_cluster_ids)}] orbit {oid}: "
            f"{n_eqs} eqs scanned; new components reached so far = "
            f"{len(new_touched)} (running)")
    log(f"  TEST 1 elapsed: {time.time()-t0:.1f} s")
    return new_touched, bridge_records


# ============================================================================
# TEST 2 -- depth-2 Laurent extensions on cluster + non-cluster orbits
# ============================================================================

def test2_depth2_bridges(rep_to_component, rep_to_id, id_to_rec,
                          cluster_ids, target_components, log,
                          time_budget_s=3600):
    """
    For each cluster + non-cluster non-tri orbit, enumerate depth-2
    fingerprints (Laurent order leading + 2) and check whether any
    triangulation cousin in those equations lies in a target_component
    (= still-untouched).

    LOGIC:  same shape as TEST 1 but at depth 2.  Iterate
            build_depth_d_fingerprints(M, zone, 9, 2), evaluate each,
            filter cousins for is_triangulation, map to components.

    PHYSICS:  depth-2 cascade equations have larger candidate pools
              and richer cousin structure than depth-1, so they may
              reach triangulation orbits that depth-1 misses.
    """
    n = 9
    all_oids = sorted(set(id_to_rec.keys()))
    log(f"TEST 2: enumerating depth-2 fingerprints for "
        f"{len(all_oids)} orbits, looking for bridges into "
        f"{len(target_components)} untouched components...")

    new_touched = set()
    bridge_records = []
    t0 = time.time()
    for i, oid in enumerate(all_oids):
        if time.time() - t0 > time_budget_s:
            log(f"  TIME BUDGET EXCEEDED at orbit {i+1}/{len(all_oids)}; "
                f"stopping TEST 2.")
            break
        M = id_to_rec[oid]["rep"]
        for zone_r in range(1, n + 1):
            cls = classify_chords_at_zone(M, zone_r, n)
            if cls["n_subs"] == 0:
                continue
            target_order = cls["ell_Z"] + 2
            fps = build_depth_d_fingerprints(M, zone_r, n, 2)
            for fp_rec in fps:
                eq = fingerprint_equation_for_M(
                    M, zone_r, n, target_order, fp_rec["fp"],
                    pool_size_cap=8000)
                if eq is None or not eq:
                    continue
                tri_cousins_in_target = []
                for ms, sc in eq:
                    if not is_triangulation(list(ms), n):
                        continue
                    rep = canonical_orbit_rep(ms, n)
                    cid = rep_to_component.get(rep)
                    if cid in target_components:
                        tri_cousins_in_target.append({
                            "tri_rep": [list(c) for c in rep],
                            "scalar": int(sc),
                            "step2_component": cid,
                        })
                        new_touched.add(cid)
                if tri_cousins_in_target:
                    bridge_records.append({
                        "from_orbit": oid,
                        "from_orbit_rep": [list(c) for c in M],
                        "zone": zone_r,
                        "depth": 2,
                        "fingerprint": fp_rec["fp_string"],
                        "tri_cousins": tri_cousins_in_target,
                    })
        if (i + 1) % 10 == 0:
            log(f"  [{i+1}/{len(all_oids)}] checked; bridges into target "
                f"so far = {len(new_touched)}; elapsed "
                f"{time.time()-t0:.1f} s")
    log(f"  TEST 2 elapsed: {time.time()-t0:.1f} s")
    return new_touched, bridge_records


# ============================================================================
# TEST 4 -- structural characterization of remaining untouched components
# ============================================================================

def cyclic_stabilizer(rep, n=9):
    """
    Smallest s > 0 such that shifting rep by s mod n leaves it
    invariant.  Stabilizer order = n / s.
    """
    from analyze_recipes import shift_multiset
    for s in range(1, n + 1):
        if shift_multiset(rep, s, n) == tuple(sorted(rep)):
            return n // s
    return 1


def chord_lengths(rep, n=9):
    """List of cyclic-distance lengths of the chords in rep."""
    out = []
    for (i, j) in rep:
        d = min(abs(j - i), n - abs(j - i))
        out.append(d)
    return sorted(out)


def ear_vertices(rep, n=9):
    """Vertices that are 'ears' of the triangulation (no chord touches)."""
    used = set()
    for (i, j) in rep:
        used.add(i)
        used.add(j)
    return sorted(set(range(1, n + 1)) - used)


def test4_structural_chars(untouched_components, rep_to_component,
                            tri_id_to_rep, log):
    """
    For each untouched component, list its triangulation orbit reps
    with structural features.

    LOGIC:  group reps by component ID; compute features.

    PHYSICS:  the goal is to identify a common feature (cyclic
              stabilizer, chord-length signature, ear-vertex pattern)
              across the still-untouched components, hinting at the
              missing mechanism.
    """
    rep_by_component = defaultdict(list)
    for rep, cid in rep_to_component.items():
        if cid in untouched_components:
            rep_by_component[cid].append(rep)

    out = {}
    for cid, reps in rep_by_component.items():
        comp_data = {"component_id": cid, "n_orbits": len(reps), "orbits": []}
        for rep in reps:
            comp_data["orbits"].append({
                "rep": [list(c) for c in rep],
                "stabilizer_order": cyclic_stabilizer(rep, 9),
                "chord_lengths": chord_lengths(rep, 9),
                "ear_vertices": ear_vertices(rep, 9),
            })
        out[cid] = comp_data
    return out


# ============================================================================
# MAIN
# ============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    log_lines = []

    def log(s):
        log_lines.append(s)
        print(s, flush=True)

    log("=" * 78)
    log("Bridge the 4 untouched Step-2 components at n=9")
    log("=" * 78)

    rep_to_id, id_to_rec = load_orbit_manifest()
    cluster_ids = load_cluster_orbits()
    log(f"\nManifest: {len(id_to_rec)} step-1 non-tri orbits at n=9")
    log(f"Cluster:  {len(cluster_ids)} orbits "
        f"(non-cluster: {len(id_to_rec) - len(cluster_ids)})")

    rep_to_component, rep_to_tri_id = load_step2_components()
    n_components = max(rep_to_component.values())
    log(f"Step-2 components on triangulations: {n_components}")
    log(f"Triangulation orbits: {len(rep_to_tri_id)}")

    initial_touched = load_initial_touched_components()
    log(f"Initial touched components (cluster ext-tri columns): "
        f"{len(initial_touched)} of {n_components}")
    untouched = set(range(1, n_components + 1)) - initial_touched
    log(f"Initial UNTOUCHED components: {sorted(untouched)}")
    log("")

    # ----- TEST 1 -----
    new_touched_1, bridges_1 = test1_non_cluster_bridges(
        rep_to_component, rep_to_id, id_to_rec, cluster_ids, log)
    initial_plus_1 = initial_touched | new_touched_1
    untouched_after_1 = set(range(1, n_components + 1)) - initial_plus_1
    log(f"\nAfter TEST 1: touched = {len(initial_plus_1)} / {n_components}; "
        f"new bridges into components {sorted(new_touched_1 - initial_touched)}")
    log(f"After TEST 1 untouched: {sorted(untouched_after_1)}")
    with open(os.path.join(OUT_DIR, "n9_bridge_test1_equations.json"),
              "w") as f:
        json.dump({
            "n_bridges": len(bridges_1),
            "new_touched_components_relative_to_initial":
                sorted(new_touched_1 - initial_touched),
            "components_touched_by_test1_total": sorted(initial_plus_1),
            "components_still_untouched": sorted(untouched_after_1),
            "bridges": bridges_1,
        }, f, indent=2)
    if not untouched_after_1:
        log("\n*** ALL 16 components covered after TEST 1.  STOP. ***")
        with open(os.path.join(OUT_DIR, "n9_bridge_test_summary.txt"),
                  "w") as f:
            f.write("\n".join(log_lines) + "\n")
        return
    log(f"  -> {len(untouched_after_1)} components remain; proceeding to TEST 2.\n")

    # ----- TEST 2 -----
    new_touched_2, bridges_2 = test2_depth2_bridges(
        rep_to_component, rep_to_id, id_to_rec, cluster_ids,
        untouched_after_1, log)
    initial_plus_2 = initial_plus_1 | new_touched_2
    untouched_after_2 = set(range(1, n_components + 1)) - initial_plus_2
    log(f"\nAfter TEST 2: touched = {len(initial_plus_2)} / {n_components}; "
        f"new bridges into components {sorted(new_touched_2)}")
    log(f"After TEST 2 untouched: {sorted(untouched_after_2)}")
    with open(os.path.join(OUT_DIR, "n9_bridge_test2_equations.json"),
              "w") as f:
        json.dump({
            "n_bridges": len(bridges_2),
            "components_touched_by_test2_relative_to_test1":
                sorted(new_touched_2),
            "components_still_untouched": sorted(untouched_after_2),
            "bridges": bridges_2,
        }, f, indent=2)
    if not untouched_after_2:
        log("\n*** ALL 16 components covered after TEST 2.  STOP. ***")
        with open(os.path.join(OUT_DIR, "n9_bridge_test_summary.txt"),
                  "w") as f:
            f.write("\n".join(log_lines) + "\n")
        return

    # ----- TEST 4 (skip TEST 3 for now -- pure-tri Laurent is heavier) -----
    log("\nTEST 4: structural characterization of remaining "
        f"{len(untouched_after_2)} untouched components.")
    structural = test4_structural_chars(
        untouched_after_2, rep_to_component, rep_to_tri_id, log)

    out_md = []
    out_md.append(f"# Structural data for {len(structural)} untouched "
                  f"Step-2 components at n=9\n")
    for cid, data in sorted(structural.items()):
        out_md.append(f"## Component {cid}  ({data['n_orbits']} orbits)\n")
        out_md.append("| Orbit rep | Stab order | Chord lengths | Ear vertices |")
        out_md.append("|---|---:|---|---|")
        for o in data["orbits"]:
            out_md.append(
                f"| `{format_multiset(tuple(tuple(c) for c in o['rep']))}` "
                f"| {o['stabilizer_order']} "
                f"| {o['chord_lengths']} "
                f"| {o['ear_vertices']} |")
        out_md.append("")
    with open(os.path.join(OUT_DIR, "n9_untouched_structural.md"),
              "w") as f:
        f.write("\n".join(out_md) + "\n")

    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)
    log(f"After TESTs 1 + 2:  {len(initial_plus_2)} / {n_components} "
        f"Step-2 components are touched.")
    log(f"Still UNTOUCHED:   {sorted(untouched_after_2)}")
    log(f"Structural characterization saved to "
        f"outputs/n9_untouched_structural.md")
    with open(os.path.join(OUT_DIR, "n9_bridge_test_summary.txt"),
              "w") as f:
        f.write("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
