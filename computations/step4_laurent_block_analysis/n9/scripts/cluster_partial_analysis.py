"""
================================================================================
cluster_partial_analysis.py  --  Fixed-cluster rank/nullity analysis at n=9
================================================================================

PURPOSE
-------
Take a FIXED set of orbit IDs as the cluster (no BFS growth), build the
depth-1 fingerprint matrix on just those orbits, and compute rank/nullity.
Any cousin that lives outside the cluster (whether it's a step-1-survivor
in a non-cluster orbit, or a triangulation, or a non-tri orbit not yet in
the manifest) becomes one column per *distinct rep* in the "external"
column block.

This lets us decide whether the discovered subsystem is already close to
full rank, without waiting for BFS to terminate.

USAGE
-----
  python3 cluster_partial_analysis.py        # uses default orbit list (90 orbits known)
  python3 cluster_partial_analysis.py 22 88 96 102 ...   # custom list

OUTPUTS
-------
  cluster_partial.json   -- machine-readable
  cluster_partial.md     -- human-readable summary

CODE STYLE: per repo standards, every function has a LOGIC + PHYSICS
docstring section.
"""

import json
import os
import sys
import time

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.normpath(os.path.join(HERE, "..", "outputs"))
DIAG_DIR = os.path.normpath(os.path.join(HERE, "..", "diagnostics"))
N8_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "n8", "scripts"))
sys.path.insert(0, N8_DIR)

from cascade_kill_n8 import (  # noqa: E402
    layer0_kill_zone, format_multiset,
)
from analyze_recipes import canonical_orbit_rep  # noqa: E402
from cluster_analysis import (  # noqa: E402
    load_orbit_manifest, depth1_fingerprint_equations,
)

# Monkey-patch: replace laurent_coefficient with a memoized version so
# repeated calls on the same (multiset, zone) pair are free.  This is
# the single biggest speedup for the matrix build, since many candidate
# multisets get evaluated multiple times across different fingerprints.
import cascade_kill_n8 as _ckn8  # noqa: E402
import diagnose_anomalies as _diag  # noqa: E402

_LAURENT_CACHE = {}
_orig_laurent_coefficient = _ckn8.laurent_coefficient


def _cached_laurent_coefficient(monomial_chords, zone_r, n, max_order):
    key = (tuple(sorted(monomial_chords)), zone_r, n, max_order)
    if key not in _LAURENT_CACHE:
        _LAURENT_CACHE[key] = _orig_laurent_coefficient(
            monomial_chords, zone_r, n, max_order)
    return _LAURENT_CACHE[key]


_ckn8.laurent_coefficient = _cached_laurent_coefficient
_diag.laurent_coefficient = _cached_laurent_coefficient


def build_partial_matrix(cluster_ids, rep_to_id, id_to_rec, n=9,
                         verbose=True):
    """
    Build the depth-1 fingerprint matrix on a fixed cluster.

    LOGIC:
      For each orbit O in cluster_ids, enumerate every K2-clean
      depth-1 fingerprint equation; for each cousin in the equation,
      add the cousin's scalar to the column corresponding to its
      orbit (cluster) or to a new "external" column (if the cousin
      lives outside the cluster).  Step-1-killed cousins contribute
      0 and drop out.

    PHYSICS:  external columns represent unknowns that the cluster's
              depth-1 equations alone cannot resolve.  These can be
              (i) non-tri step-1 survivors in orbits not in the
              cluster, (ii) triangulations (which carry tree-amplitude
              coefficients), or (iii) any other multiset whose
              canonical rep isn't in the manifest.
    """
    cluster_sorted = sorted(set(cluster_ids))
    col_index = {oid: i for i, oid in enumerate(cluster_sorted)}
    rows = []
    row_labels = []
    external_columns = {}  # canonical rep (tuple of tuples) -> col index
    external_kind = {}     # canonical rep -> "non_cluster_orbit"|"triangulation"|"unknown"

    n_cluster = len(cluster_sorted)

    def col_for_external(cousin_rep, cousin_oid_or_none):
        """Return column index for an external cousin; allocate if new."""
        if cousin_rep in external_columns:
            return external_columns[cousin_rep]
        idx = n_cluster + len(external_columns)
        external_columns[cousin_rep] = idx
        if cousin_oid_or_none is not None:
            external_kind[cousin_rep] = (
                f"non_cluster_orbit_{cousin_oid_or_none}"
            )
        else:
            external_kind[cousin_rep] = "triangulation_or_unknown"
        return idx

    for i_orbit, oid in enumerate(cluster_sorted):
        if verbose:
            print(f"  [{i_orbit+1}/{n_cluster}] orbit {oid}: enumerating "
                  f"depth-1 fingerprints...", flush=True)
        t0 = time.time()
        M = id_to_rec[oid]["rep"]
        n_eqs_this = 0
        for eq_rec in depth1_fingerprint_equations(M, n):
            row_dict = {}  # col -> scalar
            for term in eq_rec["equation"]:
                ms = tuple(tuple(c) for c in term["multiset"])
                # Step-1 killed -> contributes 0; skip.
                if layer0_kill_zone(list(ms), n) is not None:
                    continue
                cousin_rep = canonical_orbit_rep(ms, n)
                cousin_oid = rep_to_id.get(cousin_rep)
                if cousin_oid is not None and cousin_oid in col_index:
                    col = col_index[cousin_oid]
                else:
                    col = col_for_external(cousin_rep, cousin_oid)
                row_dict[col] = row_dict.get(col, 0) + term["scalar"]
            if any(v != 0 for v in row_dict.values()):
                rows.append(row_dict)
                row_labels.append({
                    "orbit": oid,
                    "zone_r": eq_rec["zone_r"],
                    "substitute": list(eq_rec["substitute"]),
                    "companion": list(eq_rec["companion"]),
                    "fingerprint": eq_rec["fingerprint"],
                })
                n_eqs_this += 1
        if verbose:
            print(f"    {n_eqs_this} eqs from orbit {oid} in "
                  f"{time.time()-t0:.1f} s; externals: {len(external_columns)}; "
                  f"cache size: {len(_LAURENT_CACHE)}", flush=True)
        # Save partial state every 5 orbits so we can inspect mid-run.
        if (i_orbit + 1) % 5 == 0 or i_orbit + 1 == n_cluster:
            n_total_cols_so_far = n_cluster + len(external_columns)
            partial = {
                "completed_orbits": i_orbit + 1,
                "matrix_so_far": [
                    [r.get(j, 0) for j in range(n_total_cols_so_far)]
                    for r in rows
                ],
                "row_labels": list(row_labels),
                "cluster_sorted": cluster_sorted,
                "external_columns": [
                    [list(c) for c in rep]
                    for rep in external_columns
                ],
                "n_cluster": n_cluster,
                "n_external": len(external_columns),
            }
            with open(os.path.join(OUT_DIR, "cluster_partial_progress.json"),
                      "w") as f:
                json.dump(partial, f)

    n_total_cols = n_cluster + len(external_columns)
    # Materialize as full matrix for sympy.
    matrix = []
    for r in rows:
        row = [r.get(i, 0) for i in range(n_total_cols)]
        matrix.append(row)

    return {
        "cluster_sorted": cluster_sorted,
        "col_index": col_index,
        "external_columns": external_columns,
        "external_kind": external_kind,
        "matrix": matrix,
        "row_labels": row_labels,
        "n_cluster": n_cluster,
        "n_external": len(external_columns),
        "n_total_cols": n_total_cols,
    }


def compute_ranks(matrix, n_cluster, n_total_cols):
    """
    Compute rank of the full matrix and rank of the cluster-only sub-matrix.

    LOGIC:
      - full rank: sympy.Matrix(matrix).rank()
      - cluster-only rank: same on matrix[:, :n_cluster] (slice cols).

    PHYSICS:
      Full rank = constraint rank with externals as free unknowns.
      Cluster-only rank = constraint rank when externals are forced to 0
      (i.e., assume their orbits eventually die independently).
    """
    if not matrix:
        return {"full_rank": 0, "cluster_rank": 0,
                "full_nullity": n_total_cols, "cluster_nullity": n_cluster,
                "cluster_nullspace": []}
    M_full = sp.Matrix(matrix)
    full_rank = M_full.rank()
    # Cluster-only sub-matrix.
    M_cluster = M_full[:, :n_cluster]
    cluster_rank = M_cluster.rank()
    cluster_nullspace = M_cluster.nullspace()
    cluster_nullspace_lists = [
        [sp.nsimplify(v[i, 0]) for i in range(n_cluster)]
        for v in cluster_nullspace
    ]
    return {
        "full_rank": int(full_rank),
        "full_nullity": int(n_total_cols - full_rank),
        "cluster_rank": int(cluster_rank),
        "cluster_nullity": int(n_cluster - cluster_rank),
        "cluster_nullspace": cluster_nullspace_lists,
    }


def main():
    # Default: use the 90 orbits already discovered by the BFS run.
    DEFAULT_CLUSTER = [
        1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 16, 17, 19, 21, 22, 23,
        24, 25, 27, 28, 34, 35, 36, 37, 38, 39, 42, 43, 44, 45, 46, 47, 49,
        50, 51, 52, 53, 54, 55, 56, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
        69, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86,
        87, 88, 89, 90, 91, 93, 95, 96, 97, 98, 100, 102, 103, 104, 107,
        108, 110, 111, 112, 113,
    ]
    if len(sys.argv) > 1:
        cluster_ids = [int(x) for x in sys.argv[1:]]
    else:
        cluster_ids = DEFAULT_CLUSTER

    rep_to_id, id_to_rec = load_orbit_manifest()
    print(f"Cluster: {len(cluster_ids)} orbits")
    print(f"Manifest: {len(id_to_rec)} orbits total")

    t0 = time.time()
    info = build_partial_matrix(cluster_ids, rep_to_id, id_to_rec, n=9,
                                 verbose=True)
    print(f"\nBuilt matrix: {len(info['matrix'])} rows, "
          f"{info['n_total_cols']} cols ({info['n_cluster']} cluster + "
          f"{info['n_external']} external) "
          f"in {time.time()-t0:.1f} s.")

    print("\nComputing ranks...")
    ranks = compute_ranks(info["matrix"], info["n_cluster"],
                          info["n_total_cols"])
    print(f"  Full rank          = {ranks['full_rank']}")
    print(f"  Full nullity       = {ranks['full_nullity']}")
    print(f"  Cluster rank       = {ranks['cluster_rank']}")
    print(f"  Cluster nullity    = {ranks['cluster_nullity']}")
    print(f"  External cols      = {info['n_external']}")
    print(f"  Cluster size m     = {info['n_cluster']}")
    print(f"  Cluster residual r = m - cluster_rank = "
          f"{info['n_cluster'] - ranks['cluster_rank']}")

    if ranks["cluster_nullspace"]:
        print(f"\nCluster nullspace basis ({len(ranks['cluster_nullspace'])} "
              f"vectors):")
        for i, v in enumerate(ranks["cluster_nullspace"][:5]):
            nz = [(j, v[j]) for j in range(len(v)) if v[j] != 0]
            print(f"  v_{i+1}: {len(nz)} nonzero entries")
            top = sorted(nz, key=lambda p: -abs(p[1]))[:8]
            for j, val in top:
                oid = info["cluster_sorted"][j]
                print(f"    orbit {oid}: coef = {val}")

    # Save outputs.
    out_json = os.path.join(OUT_DIR, "cluster_partial.json")
    out_md = os.path.join(OUT_DIR, "cluster_partial.md")
    out_data = {
        "n_cluster_orbits": info["n_cluster"],
        "n_external_columns": info["n_external"],
        "n_rows": len(info["matrix"]),
        "n_total_cols": info["n_total_cols"],
        "ranks": ranks,
        "cluster_orbits": info["cluster_sorted"],
        "external_columns": [
            {"rep": [list(c) for c in rep],
             "kind": info["external_kind"][rep],
             "col_index": info["external_columns"][rep]}
            for rep in info["external_columns"]
        ],
        # Drop full matrix from JSON to keep file small; keep row_labels.
        "row_labels": info["row_labels"],
    }
    # ranks contains sympy-typed null-space; coerce to strings.
    out_data["ranks"]["cluster_nullspace"] = [
        [str(x) for x in v] for v in ranks["cluster_nullspace"]
    ]
    with open(out_json, "w") as f:
        json.dump(out_data, f, indent=2)

    # Markdown summary.
    lines = []
    lines.append(f"# Partial-cluster rank analysis (n=9)\n")
    lines.append(f"- Cluster orbits (FIXED, not BFS-grown): "
                 f"{info['n_cluster']} orbits")
    lines.append(f"- Depth-1 fingerprint equations: {len(info['matrix'])}")
    lines.append(f"- Total columns: {info['n_total_cols']} = "
                 f"{info['n_cluster']} cluster + {info['n_external']} external")
    lines.append("")
    lines.append(f"## Rank summary\n")
    lines.append(f"| | Rank | Nullity |")
    lines.append(f"|---|---:|---:|")
    lines.append(f"| Full matrix (cluster + external cols) | "
                 f"{ranks['full_rank']} | {ranks['full_nullity']} |")
    lines.append(f"| Cluster-only sub-matrix | "
                 f"{ranks['cluster_rank']} | {ranks['cluster_nullity']} |")
    lines.append("")
    lines.append(f"**Cluster residual r = m - cluster_rank = "
                 f"{info['n_cluster'] - ranks['cluster_rank']}**")
    lines.append("")
    lines.append(f"Reading: with all external survivors set to 0 (i.e. "
                 f"assume they die independently), the depth-1 fingerprint "
                 f"matrix on the {info['n_cluster']} cluster orbits has rank "
                 f"{ranks['cluster_rank']}, leaving "
                 f"{info['n_cluster'] - ranks['cluster_rank']} orbit "
                 f"directions unconstrained at depth-1.")
    lines.append("")
    if ranks["cluster_nullspace"]:
        lines.append(f"## Cluster nullspace top entries (first {min(5, len(ranks['cluster_nullspace']))} basis vectors)\n")
        for i, v in enumerate(ranks["cluster_nullspace"][:5]):
            nz = [(j, v[j]) for j in range(len(v)) if v[j] != 0]
            top = sorted(nz, key=lambda p: -abs(p[1]))[:10]
            lines.append(f"### Basis vector v_{i+1}  ({len(nz)} nonzero entries)\n")
            lines.append("| Orbit | Coefficient |")
            lines.append("|---:|---:|")
            for j, val in top:
                oid = info["cluster_sorted"][j]
                lines.append(f"| {oid} | {val} |")
            lines.append("")
    with open(out_md, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nWrote {out_json} and {out_md}.")


if __name__ == "__main__":
    main()
