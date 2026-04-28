"""
================================================================================
build_anomaly_subsystem.py  --  Make the n=9 anomaly system human-readable
================================================================================

PURPOSE (BIG PICTURE)
---------------------
The n=9 cascade run produced
  - 4 GENUINE failures of the singleton-cascade rule: orbits {22, 46, 88, 108}
  - 15 timeout reps that didn't finish in the 600s budget

Some of those timeouts are likely artifacts (depth-1 closes if you give
them more time -- spot checks on orbits 28, 34, 35 showed this); others
might join the genuine-failures set.  This script makes the picture
*readable*: for each anomaly we print a reduced, survivor-only
fingerprint equation so a mathematician can SEE which orbits couple to
which.

OUTPUTS
-------
  outputs/n9_failure_equations_readable.txt
      For each of the 4 failures, one or more reduced fingerprint
      equations of the form
          0  =  c1 * a_{orbit i1}  +  c2 * a_{orbit i2}  +  ...
      with step-1-killable terms dropped, and orbit-IDs replacing
      multisets.

  outputs/n9_timeout_status.txt
      For each of the 15 timeouts, says whether the orbit is in the
      90-orbit cluster around orbit 22.  If yes, prints one reduced
      fingerprint equation involving it.  If no, marks as
      "unresolved timeout-only".

  outputs/n9_anomaly_subsystem_matrix.json
      A small matrix on the union {failures + timeouts + cousin orbits
      that appear in their equations}, with rank, nullity, basis.

USAGE
-----
  python3 build_anomaly_subsystem.py

This script reads existing diagnose_orbit_*.md, cluster_partial.json,
results_cascade_n9_reps.json and orbits_n9.json -- it does not run any
new sympy work.
"""

import json
import os
import re
import sys
from collections import defaultdict

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.normpath(os.path.join(HERE, "..", "outputs"))
DIAG_DIR = os.path.normpath(os.path.join(HERE, "..", "diagnostics"))
N8_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "n8", "scripts"))
sys.path.insert(0, N8_DIR)

from cascade_kill_n8 import layer0_kill_zone  # noqa: E402
from analyze_recipes import canonical_orbit_rep  # noqa: E402


def load_orbit_manifest():
    """
    Load the orbit manifest and build a rep -> orbit_id mapping.

    LOGIC:  parse orbits_n9.json; key the canonical reps as
            tuple-of-tuples for hashability.

    PHYSICS:  every n=9 step-1 non-tri survivor multiset has a unique
              canonical orbit rep under the Z_9 cyclic action; this
              mapping lets us label cousins by orbit ID rather than
              by raw multiset.
    """
    with open(os.path.join(OUT_DIR, "orbits_n9.json")) as f:
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


def load_cascade_results():
    """
    Load the per-rep cascade results and identify failures + timeouts.

    LOGIC:  parse results_cascade_n9_reps.json; bucket records by
            (recipe_found, timed_out) status.

    PHYSICS:  failures = singleton cascade exhausted no recipe.
              timeouts = budget ran out before exhaustive search
              completed; some may be artifacts.
    """
    with open(os.path.join(OUT_DIR, "results_cascade_n9_reps.json")) as f:
        d = json.load(f)
    failures = []
    timeouts = []
    for r in d["records"]:
        if r.get("timed_out"):
            timeouts.append(r["orbit_id"])
        elif not r.get("recipe_found"):
            failures.append(r["orbit_id"])
    return sorted(failures), sorted(timeouts)


def parse_diagnose_md(orbit_id):
    """
    Parse a diagnose_orbit_<id>.md file to extract the depth-1
    fingerprint equations of M_rep.

    LOGIC:  regex-scan the file for blocks of the form

      - Z_{r,r+2}, fp `<expr>`: M scalar=<int>, K survivor-cousins out of N
        - survivor cousin: <list> (scalar S)
        - survivor cousin: <list> (scalar S)
        ...

      Each match yields a {zone, fingerprint, M_scalar, survivors} dict.
      Step-1-killable cousins are *not* listed in the diagnose .md
      output (the script already filtered them), so this gives the
      "reduced" survivor-only equation directly.

    PHYSICS:  these are exactly the equations the user wants -- with
              step-1-killable terms already eliminated, only orbit-to-
              orbit couplings remain.
    """
    path = os.path.join(DIAG_DIR, f"diagnose_orbit_{orbit_id}.md")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        text = f.read()
    out = []
    # Find blocks starting with "- Z_{a,b}, fp `...`: M scalar=...,
    # K survivor-cousins out of N"
    pattern = re.compile(
        r"-\s*Z_\{(\d+),(\d+)\},\s*fp\s*`([^`]+)`:\s*M scalar=(-?\d+),"
        r"\s*(\d+) survivor-cousins out of (\d+)\n"
        r"((?:  -\s*survivor cousin:.+\n)*)"
    )
    for m in pattern.finditer(text):
        zone = (int(m.group(1)), int(m.group(2)))
        fp = m.group(3)
        m_scalar = int(m.group(4))
        n_survivors = int(m.group(5))
        n_total_cousins = int(m.group(6))
        survivor_block = m.group(7)
        survivors = []
        for ml in re.finditer(
            r"survivor cousin:\s*(\[\[.+?\]\])\s*\(scalar (-?\d+)\)",
            survivor_block,
        ):
            ms_str = ml.group(1)
            scalar = int(ml.group(2))
            ms = json.loads(ms_str)
            survivors.append({
                "multiset": [tuple(c) for c in ms],
                "scalar": scalar,
            })
        out.append({
            "zone": zone,
            "fingerprint": fp,
            "M_scalar": m_scalar,
            "n_survivors": n_survivors,
            "n_total_cousins": n_total_cousins,
            "survivors": survivors,
        })
    return out


def format_orbit_eq(orbit_id, M_rep, eq, rep_to_id, id_to_rec):
    """
    Format one fingerprint equation as a human-readable line:
        0  =  c0 a_{O22} + c1 a_{O88} - c2 a_{O96} + ...

    LOGIC:  the equation has M (= the orbit being analysed) plus
            survivor-cousin orbits.  Resolve each survivor multiset
            to its canonical orbit ID via rep_to_id.  Render terms
            in sorted orbit-ID order.

    PHYSICS:  by cyclic symmetry, all multisets in the same orbit
              share a coefficient, so listing by orbit ID is sufficient.
              Step-1-killable terms have already been removed in the
              diagnose .md, so this IS the reduced equation.
    """
    terms = []
    # The orbit M itself
    terms.append((orbit_id, eq["M_scalar"]))
    for s in eq["survivors"]:
        ms = tuple(s["multiset"])
        rep = canonical_orbit_rep(ms, 9)
        oid = rep_to_id.get(rep)
        if oid is None:
            # Triangulation cousin or unknown: keep as raw multiset.
            terms.append((f"raw{tuple(ms)}", s["scalar"]))
        else:
            terms.append((oid, s["scalar"]))
    # Aggregate by orbit (in case multiple cousins land in same orbit).
    agg = defaultdict(int)
    for label, sc in terms:
        agg[label] += sc
    # Format.
    parts = []
    for label, sc in sorted(agg.items(),
                             key=lambda x: (isinstance(x[0], str), x[0])):
        if sc == 0:
            continue
        sign = "+" if sc > 0 else "-"
        magnitude = abs(sc)
        if magnitude == 1:
            term = f"a_{{O{label}}}" if isinstance(label, int) else f"a_{label}"
        else:
            term = f"{magnitude} a_{{O{label}}}" if isinstance(label, int) else f"{magnitude} a_{label}"
        if not parts:
            parts.append(f"{'-' if sc < 0 else ''}{term}")
        else:
            parts.append(f" {sign} {term}")
    return "0 = " + "".join(parts)


def cluster_orbits_from_partial():
    """
    Read the 90-orbit cluster from cluster_partial.json.

    LOGIC:  the JSON has cluster_orbits = the sorted list of orbit IDs
            in the BFS-discovered cluster.

    PHYSICS:  these are the orbits the depth-1 cousin graph BFS pulled
              in from orbit 22 before being terminated for the rank
              analysis.  They cover the cluster the failures live in.
    """
    path = os.path.join(OUT_DIR, "cluster_partial.json")
    if not os.path.exists(path):
        return set()
    with open(path) as f:
        d = json.load(f)
    return set(d.get("cluster_orbits", []))


def find_eq_involving_orbit(target_orbit_id, rep_to_id, id_to_rec,
                             diagnose_orbits=(22, 25, 28, 34, 35,
                                              46, 88, 108)):
    """
    Search the parsed diagnose .md files for any fingerprint equation
    whose survivor cousins include the target orbit.

    LOGIC:  iterate diagnose_orbits, parse each, and check if any
            survivor cousin in any equation maps to target_orbit_id.

    PHYSICS:  gives ONE example equation showing the target orbit
              coupled to its reaction.  More fingerprint equations
              involving any orbit can be built by running the cascade
              search on its rep, but for the readability output one
              example suffices.
    """
    for src_oid in diagnose_orbits:
        eqs = parse_diagnose_md(src_oid)
        if eqs is None:
            continue
        if src_oid not in id_to_rec:
            continue
        src_rep = id_to_rec[src_oid]["rep"]
        for eq in eqs:
            for s in eq["survivors"]:
                ms = tuple(s["multiset"])
                rep = canonical_orbit_rep(ms, 9)
                oid = rep_to_id.get(rep)
                if oid == target_orbit_id:
                    return src_oid, src_rep, eq
    return None


def build_subsystem_matrix(failure_ids, timeout_ids, rep_to_id, id_to_rec,
                           diagnose_orbits=(22, 25, 28, 34, 35, 46, 88, 108)):
    """
    Build a small fingerprint matrix on {failures + timeouts + cousin
    orbits appearing in their equations}, then compute rank/nullity.

    LOGIC:  start with anomaly_orbits = failure_ids + timeout_ids.
            Parse each available diagnose_orbit_<id>.md to collect rows;
            each row is one fingerprint equation, columns indexed by
            anomaly_orbits + any cousin orbits not already there.
            Compute sympy rank + nullspace.

    PHYSICS:  this is a *small* version of the cluster-partial matrix
              -- just enough orbits to cover the anomalies and their
              direct cousins.  If the small system has zero nullity,
              the anomalies are jointly killed by their fingerprint
              equations alone.
    """
    cols = list(failure_ids) + [t for t in timeout_ids
                                 if t not in failure_ids]
    cols_set = set(cols)
    rows = []
    row_labels = []
    extra_orbit_cols = []
    for src_oid in diagnose_orbits:
        eqs = parse_diagnose_md(src_oid)
        if not eqs:
            continue
        for eq in eqs:
            row_dict = defaultdict(int)
            row_dict[src_oid] += eq["M_scalar"]
            for s in eq["survivors"]:
                ms = tuple(s["multiset"])
                rep = canonical_orbit_rep(ms, 9)
                oid = rep_to_id.get(rep)
                if oid is None:
                    # External (triangulation or unknown).  Use a
                    # synthetic column id (negative integer to avoid
                    # collision).
                    key = -hash(rep) & 0xFFFFFFFF
                    row_dict[key] += s["scalar"]
                else:
                    row_dict[oid] += s["scalar"]
            # Extend columns as new orbit-ids appear.
            for k in list(row_dict.keys()):
                if k not in cols_set:
                    cols.append(k)
                    cols_set.add(k)
                    extra_orbit_cols.append(k)
            row = [row_dict[c] for c in cols]
            if any(v != 0 for v in row):
                rows.append(row)
                row_labels.append({
                    "from_orbit": src_oid,
                    "zone": list(eq["zone"]),
                    "fingerprint": eq["fingerprint"],
                })
    n_total_cols = len(cols)
    # Pad rows to current column count (since we extended cols mid-loop).
    rows = [r + [0] * (n_total_cols - len(r)) for r in rows]

    if not rows:
        return {"cols": cols, "rows": rows, "rank": 0, "nullity": n_total_cols,
                "row_labels": row_labels}

    Mat = sp.Matrix(rows)
    rank = int(Mat.rank())
    nullity = n_total_cols - rank
    null_basis = [
        [sp.nsimplify(v[i, 0]) for i in range(n_total_cols)]
        for v in Mat.nullspace()
    ]
    return {
        "cols": cols,
        "rows": rows,
        "rank": rank,
        "nullity": nullity,
        "null_basis": null_basis,
        "row_labels": row_labels,
        "extra_orbit_cols": extra_orbit_cols,
    }


def main():
    rep_to_id, id_to_rec = load_orbit_manifest()
    failures, timeouts = load_cascade_results()
    cluster = cluster_orbits_from_partial()
    print(f"Failures: {failures}")
    print(f"Timeouts: {timeouts}")
    print(f"Cluster (from cluster_partial.json): {len(cluster)} orbits")
    print()

    # ------------------------------------------------------------------
    # 1.  Failure equations (readable)
    # ------------------------------------------------------------------
    out_path = os.path.join(OUT_DIR, "n9_failure_equations_readable.txt")
    with open(out_path, "w") as f:
        f.write("=" * 78 + "\n")
        f.write("Reduced (survivor-only) depth-1 fingerprint equations\n"
                "for the 4 genuine n=9 cascade-failure orbits.\n")
        f.write("Step-1-killable cousins have been dropped (they "
                "contribute 0).\n")
        f.write("=" * 78 + "\n\n")
        for oid in failures:
            f.write("-" * 78 + "\n")
            rec = id_to_rec.get(oid)
            if rec is None:
                f.write(f"Orbit {oid}: not in manifest.\n\n")
                continue
            f.write(f"Orbit {oid}  M_rep = "
                    f"{rec['rep']}  V_missed = "
                    f"{rec['vertices_missed']}\n")
            f.write("-" * 78 + "\n")
            eqs = parse_diagnose_md(oid)
            if not eqs:
                f.write(f"  (no diagnose_orbit_{oid}.md found)\n\n")
                continue
            for eq in eqs:
                f.write(f"\n  Z_{{{eq['zone'][0]},{eq['zone'][1]}}}  "
                        f"fp = {eq['fingerprint']}\n")
                f.write(f"    {format_orbit_eq(oid, rec['rep'], eq, rep_to_id, id_to_rec)}\n")
                f.write(f"    ({eq['n_survivors']} survivor-cousins out of "
                        f"{eq['n_total_cousins']} total cousins; "
                        f"the rest are step-1-killable.)\n")
            f.write("\n")
    print(f"Wrote {out_path}")

    # ------------------------------------------------------------------
    # 2.  Timeout status
    # ------------------------------------------------------------------
    out_path = os.path.join(OUT_DIR, "n9_timeout_status.txt")
    with open(out_path, "w") as f:
        f.write("=" * 78 + "\n")
        f.write("Status of the 15 n=9 timeout-orbit reps.\n")
        f.write("Cluster membership refers to the 90-orbit BFS cluster\n"
                "around orbit 22 (cluster_partial.json).\n")
        f.write("=" * 78 + "\n\n")
        for oid in timeouts:
            rec = id_to_rec.get(oid, {"rep": "?", "vertices_missed": "?"})
            in_cluster = oid in cluster
            f.write("-" * 78 + "\n")
            f.write(f"Orbit {oid}  M_rep = {rec.get('rep')}  "
                    f"V_missed = {rec.get('vertices_missed')}\n")
            f.write(f"  In 90-orbit cluster? {'YES' if in_cluster else 'NO'}\n")
            if in_cluster:
                hit = find_eq_involving_orbit(oid, rep_to_id, id_to_rec)
                if hit:
                    src_oid, src_rep, eq = hit
                    f.write(f"  Example equation (from orbit {src_oid}'s "
                            f"depth-1 search):\n")
                    f.write(f"    Z_{{{eq['zone'][0]},{eq['zone'][1]}}}  "
                            f"fp = {eq['fingerprint']}\n")
                    f.write(f"    {format_orbit_eq(src_oid, src_rep, eq, rep_to_id, id_to_rec)}\n")
                else:
                    f.write(f"  No example equation cached "
                            f"(would need to run cascade search on "
                            f"orbit {oid}'s rep).\n")
            else:
                f.write(f"  UNRESOLVED-TIMEOUT-ONLY: cluster expansion "
                        f"would need more BFS to confirm recipe.\n")
            f.write("\n")
    print(f"Wrote {out_path}")

    # ------------------------------------------------------------------
    # 3.  Anomaly subsystem matrix
    # ------------------------------------------------------------------
    info = build_subsystem_matrix(failures, timeouts, rep_to_id, id_to_rec)
    out_path = os.path.join(OUT_DIR, "n9_anomaly_subsystem_matrix.json")
    serial = {
        "cols_orbit_ids_or_external_hashes": info["cols"],
        "n_cols": len(info["cols"]),
        "n_rows": len(info["rows"]),
        "rank": info["rank"],
        "nullity": info["nullity"],
        "row_labels": info["row_labels"],
        "extra_orbit_cols": info.get("extra_orbit_cols", []),
        "null_basis": [[str(x) for x in v]
                        for v in info.get("null_basis", [])],
    }
    with open(out_path, "w") as f:
        json.dump(serial, f, indent=2)
    print(f"Wrote {out_path}")
    print(f"  matrix dims: {len(info['rows'])} rows x "
          f"{len(info['cols'])} cols, rank {info['rank']}, "
          f"nullity {info['nullity']}")


if __name__ == "__main__":
    main()
