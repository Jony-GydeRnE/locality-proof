"""
================================================================================
cluster_analysis.py  --  Coupled-cascade cluster analysis at n=9
================================================================================

PURPOSE (BIG PICTURE)
---------------------
At n=9 the depth-1 Laurent-cascade kill mechanism that closed every orbit at
n=8 fails on orbit 22 -- not because of a search bug or a missing depth, but
because the cousins of M_{22} in its depth-1 fingerprint equations are
themselves step-1 survivors.  This is the "coupled subsystem" scenario:
orbit 22 cannot be killed in isolation; it must be killed jointly with the
other orbits its cousins live in.

This script does the cluster-and-rank analysis the user asked for:

  Step 1.  Cluster identification.  BFS from orbit 22: an orbit B is added
           to the cluster if some orbit-A member's depth-1 fingerprint
           equation has B's representative as a step-1-survivor cousin.
           Iterate until closed.

  Step 2.  Cluster matrix M_FP.  For each cluster orbit O and each
           (zone Z, substitute U) with K2-clean U, build the depth-1
           fingerprint equation; map each cousin to the canonical orbit
           rep; sum scalars per orbit; record one row of M_FP whose
           columns index the cluster orbits.  Step-1-killed cousins
           contribute 0 (they're known zero) so they drop out.  Cousins
           in non-cluster step-1-survivor orbits should not appear if
           the cluster is closed (sanity check).

  Step 3.  Rank & nullspace.  Compute over Q via sympy.  Output rank,
           residual r = m - rank, and a basis for the nullspace.

  Step 4.  Anchor search.  For each residual basis direction, search
           depth-2 (and depth-3 if needed) fingerprints across all
           zones of all cluster orbits, looking for an equation that
           kills the remaining direction (after applying the M_FP
           reductions).

  Step 5.  Compare to the n=6 perfect-matching cluster.  At n=6, the
           4 perfect matchings form a cluster with M_FP rank 3 and a
           single anchor.  We report whether C_22 mirrors that
           structure or is genuinely larger.

OUTPUTS
-------
  cluster_22.json            -- cluster membership, M_FP rows, rank,
                                 nullspace basis, anchor recipes (if any).
  cluster_22_summary.md      -- human-readable summary.
"""

import json
import os
import sys
import time
from collections import defaultdict, deque
from itertools import combinations_with_replacement

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
from analyze_recipes import (  # noqa: E402
    canonical_orbit_rep, vertices_of, missed_vertices,
)
from diagnose_anomalies import (  # noqa: E402
    classify_chords_at_zone, build_depth_d_fingerprints,
    fingerprint_equation_for_M, free_vars_at_zone,
)


# Globals: orbit manifest helpers.
N = 9


def load_orbit_manifest(path=None):
    """Load orbits_n9.json and return rep_to_id, id_to_record."""
    if path is None:
        path = os.path.join(OUT_DIR, "orbits_n9.json")
    with open(path) as f:
        d = json.load(f)
    rep_to_id = {}
    id_to_rec = {}
    for o in d["orbits"]:
        rep = tuple(tuple(c) for c in o["representative"])
        rep_to_id[rep] = o["orbit_id"]
        id_to_rec[o["orbit_id"]] = {
            "rep": rep,
            "size": o["size"],
            "vertices_missed": o["vertices_missed"],
        }
    return rep_to_id, id_to_rec


# ============================================================================
# STEP 1 -- CLUSTER IDENTIFICATION
# ============================================================================

def depth1_fingerprint_equations(M, n):
    """
    Yield (zone_r, substitute U, companion Y_U, fingerprint, equation) for
    every K2-clean (Z, U) at depth-1 (Laurent order leading + 1) where
    M's representative is in the equation with nonzero scalar.
    """
    M_set = set(M)
    M_sorted = tuple(sorted(M))
    for zone_r in range(1, n + 1):
        cls = classify_chords_at_zone(M, zone_r, n)
        if cls["n_subs"] == 0:
            continue
        for c, kind, comp in cls["bucket"]:
            if kind != "substitute":
                continue
            if comp in M_set:
                continue  # K2 violation -- skip this (Z, U)
            U = c
            Y = comp
            # Build depth-1 fingerprint:  X_Y / (rest free chord factors).
            ms_minus = list(M)
            ms_minus.remove(U)
            sub_set_local = {sub for _, sub in cls["pairs"]}
            rest_free = [
                cc for cc in ms_minus
                if cc != cls["special"] and cc != cls["bare"]
                and cc not in sub_set_local
            ]
            fp = chord_symbol(Y)
            for cc in rest_free:
                fp /= chord_symbol(cc)
            target_order = cls["ell_Z"] + 1
            eq = fingerprint_equation_for_M(M, zone_r, n, target_order,
                                             sp.simplify(fp),
                                             pool_size_cap=10000)
            if eq is None or not eq:
                continue
            fish_scalar = None
            for m, s in eq:
                if tuple(sorted(m)) == M_sorted:
                    fish_scalar = s
                    break
            if fish_scalar is None or fish_scalar == 0:
                continue
            yield {
                "zone_r": zone_r,
                "substitute": U,
                "companion": Y,
                "fingerprint": str(sp.simplify(fp)),
                "fp_expr": fp,
                "target_order": target_order,
                "M_scalar": int(fish_scalar),
                "equation": [(tuple(c) for c in (m, s)) for m, s in eq] and [
                    {"multiset": [list(c) for c in m], "scalar": int(s)}
                    for m, s in eq
                ],
            }


def cluster_bfs(start_orbit_id, rep_to_id, id_to_rec, n=9, verbose=True):
    """
    BFS from start_orbit_id: an orbit B is added to the cluster if some
    orbit-A member's depth-1 fingerprint equation contains B's rep (or any
    member of B) as a step-1-survivor cousin.

    Returns: cluster (set of orbit_ids), edges (list of (A, B, recipe)).
    """
    cluster = {start_orbit_id}
    queue = deque([start_orbit_id])
    edges = []
    eqs_per_orbit = {}  # orbit_id -> list of fingerprint records

    while queue:
        oid = queue.popleft()
        if verbose:
            print(f"  BFS visiting orbit {oid} (rep "
                  f"{format_multiset(id_to_rec[oid]['rep'])})...", flush=True)
        M = id_to_rec[oid]["rep"]
        my_eqs = list(depth1_fingerprint_equations(M, n))
        eqs_per_orbit[oid] = my_eqs
        for eq_rec in my_eqs:
            for term in eq_rec["equation"]:
                ms = tuple(tuple(c) for c in term["multiset"])
                if tuple(sorted(ms)) == tuple(sorted(M)):
                    continue  # the orbit's own member
                if layer0_kill_zone(list(ms), n) is not None:
                    continue  # step-1 killed
                # Otherwise this cousin is a step-1 survivor; find its orbit.
                cousin_rep = canonical_orbit_rep(ms, n)
                cousin_oid = rep_to_id.get(cousin_rep)
                if cousin_oid is None:
                    if verbose:
                        print(f"    WARNING: cousin {ms} not in orbit "
                              f"manifest (canonical {cousin_rep}).")
                    continue
                edges.append({
                    "from_orbit": oid,
                    "to_orbit": cousin_oid,
                    "zone_r": eq_rec["zone_r"],
                    "fingerprint": eq_rec["fingerprint"],
                    "scalar": term["scalar"],
                })
                if cousin_oid not in cluster:
                    cluster.add(cousin_oid)
                    queue.append(cousin_oid)
                    if verbose:
                        print(f"    -> added orbit {cousin_oid} to cluster.")
    return cluster, edges, eqs_per_orbit


# ============================================================================
# STEP 2 -- BUILD M_FP MATRIX
# ============================================================================

def build_cluster_matrix(cluster, rep_to_id, id_to_rec, eqs_per_orbit, n=9,
                         verbose=True):
    """
    Build the M_FP matrix:
        rows = depth-1 fingerprint equations (one per (orbit, Z, U)),
        cols = cluster orbit ids (in sorted order).
    Entry M_FP[i, j] = sum over cousins of equation i in orbit j of the
                        cousin scalar.

    Step-1-killed cousins are ignored (their column would be a "killed"
    column with all entries 0 anyway -- but we don't include such
    columns).

    Cousins in non-cluster step-1-SURVIVOR orbits should not occur if the
    cluster is closed; we sanity-check and warn if they do.
    """
    cluster_sorted = sorted(cluster)
    col_index = {oid: i for i, oid in enumerate(cluster_sorted)}
    rows = []
    row_labels = []
    leak_warnings = []

    for oid in cluster_sorted:
        for eq_rec in eqs_per_orbit.get(oid, []):
            row = [0] * len(cluster_sorted)
            for term in eq_rec["equation"]:
                ms = tuple(tuple(c) for c in term["multiset"])
                if layer0_kill_zone(list(ms), n) is not None:
                    continue  # killed -> 0
                cousin_rep = canonical_orbit_rep(ms, n)
                cousin_oid = rep_to_id.get(cousin_rep)
                if cousin_oid is None:
                    leak_warnings.append({
                        "from_orbit": oid,
                        "multiset": [list(c) for c in ms],
                    })
                    continue
                if cousin_oid not in col_index:
                    leak_warnings.append({
                        "from_orbit": oid,
                        "to_orbit_outside_cluster": cousin_oid,
                        "multiset": [list(c) for c in ms],
                    })
                    continue
                row[col_index[cousin_oid]] += term["scalar"]
            if any(v != 0 for v in row):
                rows.append(row)
                row_labels.append({
                    "orbit": oid,
                    "zone_r": eq_rec["zone_r"],
                    "substitute": list(eq_rec["substitute"]),
                    "companion": list(eq_rec["companion"]),
                    "fingerprint": eq_rec["fingerprint"],
                })
    return rows, row_labels, cluster_sorted, leak_warnings


def matrix_rank_nullspace(rows, n_cols):
    """
    Compute rank and a nullspace basis of the integer matrix `rows` (m x
    n_cols).  Uses sympy for exact-rational arithmetic.
    """
    if not rows:
        return 0, [list(range(n_cols))]  # whole space is null
    M = sp.Matrix(rows)
    rank = M.rank()
    null_basis = M.nullspace()
    # Convert sympy column vectors to plain Python lists.
    null_basis_list = []
    for v in null_basis:
        null_basis_list.append([sp.nsimplify(v[i, 0]) for i in range(n_cols)])
    return rank, null_basis_list


# ============================================================================
# STEP 3 -- ANCHOR SEARCH (depth 2+)
# ============================================================================

def search_anchor_for_residual(cluster_sorted, id_to_rec, basis_vec, n=9,
                               max_depth=2, verbose=True):
    """
    Try to find a depth-d (d >= 2) fingerprint equation across any cluster
    orbit's zones whose non-cluster cousins are all step-1-killable AND
    whose contribution to the residual basis vector is nonzero.

    LOGIC:  for each cluster orbit O, for each zone Z, for each
            depth-d fingerprint of M_O at Z, evaluate the equation;
            project onto cluster cousins; check non-cluster cousins
            killability; check that the projected row has nonzero dot
            product with `basis_vec` (so it actually constrains the
            residual direction).

    Returns the first such anchor recipe (or None).
    """
    col_index = {oid: i for i, oid in enumerate(cluster_sorted)}
    for oid in cluster_sorted:
        M = id_to_rec[oid]["rep"]
        for depth in range(2, max_depth + 1):
            fps_per_zone = []
            for zone_r in range(1, n + 1):
                fps = build_depth_d_fingerprints(M, zone_r, n, depth)
                for fp_rec in fps:
                    fps_per_zone.append((zone_r, fp_rec))
            for zone_r, fp_rec in fps_per_zone:
                target_order = classify_chords_at_zone(M, zone_r, n)["ell_Z"] + depth
                eq = fingerprint_equation_for_M(
                    M, zone_r, n, target_order, fp_rec["fp"],
                    pool_size_cap=8000)
                if eq is None or not eq:
                    continue
                M_sorted = tuple(sorted(M))
                fish_scalar = None
                for m, s in eq:
                    if tuple(sorted(m)) == M_sorted:
                        fish_scalar = s
                        break
                if fish_scalar is None or fish_scalar == 0:
                    continue
                # Build the projected row.
                row = [0] * len(cluster_sorted)
                non_cluster_survivor = False
                for m, s in eq:
                    if layer0_kill_zone(list(m), n) is not None:
                        continue
                    cousin_rep = canonical_orbit_rep(m, n)
                    cousin_oid = None
                    for cid, rec in id_to_rec.items():
                        if rec["rep"] == cousin_rep:
                            cousin_oid = cid
                            break
                    if cousin_oid is None or cousin_oid not in col_index:
                        non_cluster_survivor = True
                        break
                    row[col_index[cousin_oid]] += int(s)
                if non_cluster_survivor:
                    continue
                # Check that this row has nonzero dot product with the
                # residual basis vec.
                dot = sum(row[i] * basis_vec[i] for i in range(len(row)))
                if dot != 0:
                    return {
                        "orbit": oid,
                        "depth": depth,
                        "zone_r": zone_r,
                        "fingerprint": fp_rec["fp_string"],
                        "row": row,
                        "dot_with_residual": int(dot),
                    }
    return None


# ============================================================================
# STEP 4 -- DRIVER
# ============================================================================

def main(start_orbit_id=22, output_json=None, output_md=None):
    rep_to_id, id_to_rec = load_orbit_manifest()
    print(f"Manifest: {len(id_to_rec)} orbits.")
    print(f"Starting cluster BFS from orbit {start_orbit_id}...")
    t0 = time.time()
    cluster, edges, eqs_per_orbit = cluster_bfs(
        start_orbit_id, rep_to_id, id_to_rec, n=N, verbose=True)
    print(f"BFS done in {time.time()-t0:.1f} s.  Cluster size: {len(cluster)}")
    print(f"Cluster orbits: {sorted(cluster)}")

    # Step 2: Matrix.
    print("\nBuilding M_FP matrix...")
    rows, row_labels, cluster_sorted, leak_warnings = build_cluster_matrix(
        cluster, rep_to_id, id_to_rec, eqs_per_orbit, n=N, verbose=True)
    print(f"  rows: {len(rows)},  cols: {len(cluster_sorted)}")
    if leak_warnings:
        print(f"  WARNING: {len(leak_warnings)} cluster-leak events "
              f"(cousins outside cluster but step-1 survivors).")

    # Step 3: Rank.
    print("\nComputing rank...")
    rank, null_basis = matrix_rank_nullspace(rows, len(cluster_sorted))
    residual = len(cluster_sorted) - rank
    print(f"  rank(M_FP) = {rank}")
    print(f"  cluster size m = {len(cluster_sorted)}")
    print(f"  residual r = m - rank = {residual}")
    print(f"  nullspace basis vectors:")
    for v in null_basis:
        print(f"    {v}")

    # Step 4: Anchors.
    anchors = []
    if residual > 0:
        print("\nSearching for anchors at depth 2 ...")
        for i, basis_vec in enumerate(null_basis):
            print(f"  Residual direction {i+1}/{len(null_basis)}: {basis_vec}")
            t1 = time.time()
            anc = search_anchor_for_residual(
                cluster_sorted, id_to_rec, basis_vec, n=N,
                max_depth=2, verbose=True)
            print(f"  ... ({time.time()-t1:.1f} s)")
            anchors.append(anc)
            if anc:
                print(f"    FOUND anchor: orbit {anc['orbit']}, depth "
                      f"{anc['depth']}, zone Z_{{{anc['zone_r']},*}}, "
                      f"fp `{anc['fingerprint']}` "
                      f"(dot with residual = {anc['dot_with_residual']})")
            else:
                print(f"    NO depth-2 anchor found for this direction.")
    else:
        print("Rank equals cluster size; no residual; cluster fully "
              "killed by depth-1 alone.")

    # Output.
    out = {
        "start_orbit": start_orbit_id,
        "cluster": cluster_sorted,
        "cluster_size": len(cluster_sorted),
        "edges_count": len(edges),
        "M_FP_rows": rows,
        "M_FP_row_labels": row_labels,
        "M_FP_columns": cluster_sorted,
        "rank": rank,
        "residual": residual,
        "nullspace_basis": [
            [str(x) for x in v] for v in null_basis
        ],
        "leak_warnings": leak_warnings,
        "anchors": anchors,
    }
    if output_json is None:
        output_json = os.path.join(OUT_DIR, f"cluster_{start_orbit_id}.json")
    if output_md is None:
        output_md = os.path.join(HERE, "..", f"cluster_{start_orbit_id}_summary.md")
    with open(output_json, "w") as f:
        json.dump(out, f, indent=2, default=str)

    # Markdown summary.
    lines = []
    lines.append(f"# Cluster analysis starting from orbit {start_orbit_id}\n")
    lines.append(f"Cluster size m = **{len(cluster_sorted)}** orbits: "
                 f"`{cluster_sorted}`.\n")
    for oid in cluster_sorted:
        rec = id_to_rec[oid]
        lines.append(f"- Orbit {oid} (size {rec['size']}): "
                     f"`{format_multiset(rec['rep'])}` "
                     f"V_missed={rec['vertices_missed']}")
    lines.append("")
    lines.append(f"## M_FP matrix\n")
    lines.append(f"- rows: {len(rows)} depth-1 fingerprint equations")
    lines.append(f"- cols: {len(cluster_sorted)} cluster orbits")
    lines.append(f"- rank(M_FP) = **{rank}**")
    lines.append(f"- residual r = m - rank = **{residual}**")
    if residual > 0:
        lines.append("\n### Nullspace basis\n")
        for i, v in enumerate(null_basis):
            lines.append(f"- v_{i+1}: " + ", ".join(
                f"orbit {cluster_sorted[j]}: {v[j]}" for j in range(len(v))
                if v[j] != 0
            ))
    if leak_warnings:
        lines.append(f"\n### Cluster-leak warnings ({len(leak_warnings)})\n")
        for w in leak_warnings[:10]:
            lines.append(f"- {w}")
    lines.append(f"\n## Anchors\n")
    if not anchors:
        lines.append("- No residual; no anchor needed.")
    for i, a in enumerate(anchors):
        if a is None:
            lines.append(f"- Direction {i+1}: NO depth-2 anchor found.")
        else:
            lines.append(f"- Direction {i+1}: anchor at orbit {a['orbit']}, "
                         f"depth {a['depth']}, zone Z_{{{a['zone_r']},*}}, "
                         f"fp `{a['fingerprint']}`, dot = "
                         f"{a['dot_with_residual']}.")

    with open(output_md, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"\nWrote {output_json} and {output_md}.")
    return out


if __name__ == "__main__":
    start_id = 22
    if len(sys.argv) > 1:
        start_id = int(sys.argv[1])
    main(start_orbit_id=start_id)
