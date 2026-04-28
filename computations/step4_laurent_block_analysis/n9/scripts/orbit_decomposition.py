"""
================================================================================
orbit_decomposition.py  --  Cyclic-orbit decomposition of n=9 step-1 survivors
================================================================================

PURPOSE (BIG PICTURE)
---------------------
At n=9 the layer-0 ("Step-1") kill mechanism leaves 1 011 non-triangulation
multisets uncaught.  Running the full Laurent cascade on all 1 011 is
expensive (~3 hours of sympy work).  But the kill mechanism is cyclically
equivariant: if M dies via a depth-1 cascade with recipe (Z, k, fp), then
its cyclic shift M' = M + s dies via the shifted recipe (Z + s, k, fp + s).
So we only need to verify the cascade once per Z_9-orbit; every other
orbit member inherits.

This script:

  (1) Enumerates the 1 011 step-1 non-triangulation survivors at n=9.
  (2) Partitions them into cyclic orbits under the Z_9 action
      v -> v + s (mod 9).
  (3) Outputs one canonical representative per orbit, plus orbit size
      (which must divide 9, so size in {1, 3, 9}).
  (4) Writes a markdown summary and a JSON manifest of orbit
      representatives for the cascade-on-reps script to consume.

It does NOT attempt any Laurent-cascade verification -- that's a separate
step in `cascade_kill_orbit_reps_n9.py`.

OUTPUTS
-------
  orbits_n9.json    -- list of {orbit_id, size, representative,
                                  vertices_used, vertices_missed,
                                  one_member_index}
  orbits_n9.md      -- human-readable summary table.

CODE STYLE
----------
Per repo contribution standards: every function has a LOGIC docstring
section and a PHYSICS section.  Most utilities are imported from
cascade_kill_n8.py (one folder up) since they are parametric in n.
"""

import json
import os
import sys
from collections import defaultdict

# Pull machinery from the cascade_n8 folder.  These functions are all
# parametric in n.
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.normpath(os.path.join(HERE, "..", "outputs"))
DIAG_DIR = os.path.normpath(os.path.join(HERE, "..", "diagnostics"))
N8_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "n8", "scripts"))
sys.path.insert(0, N8_DIR)

from cascade_kill_n8 import (  # noqa: E402
    all_chords, layer0_kill_zone, is_triangulation,
    all_size_N_multisets, format_multiset,
)
from analyze_recipes import (  # noqa: E402
    vertices_of, missed_vertices, canonical_orbit_rep,
    shift_multiset,
)


# ============================================================================
# 1.  Enumerate n=9 step-1 survivors
# ============================================================================

def step1_survivors(n=9):
    """
    Compute non-triangulation multisets of size n-3 that survive Step-1
    at every cyclic zone.

    LOGIC:  iterate every size-(n-3) multiset of chords; keep those for
            which layer0_kill_zone returns None and that fail the
            is_triangulation test.

    PHYSICS:  these are the "fish" at n -- the input set the Laurent
              cascade must handle.  At n = 9 there are 1 011 of them
              (cf. headline numbers in computations/step1_layer0_kill/).
    """
    out = []
    for ms in all_size_N_multisets(n):
        if layer0_kill_zone(list(ms), n) is None:
            if not is_triangulation(list(ms), n):
                out.append(ms)
    return out


# ============================================================================
# 2.  Orbit decomposition
# ============================================================================

def orbit_decomposition(survivors, n=9):
    """
    Group survivors into Z_n-orbits keyed by canonical representative.

    LOGIC:  for each survivor M, compute canonical_orbit_rep(M) (the
            lexicographically smallest cyclic shift), then group.

    PHYSICS:  two survivors are in the same orbit iff they differ by a
              cyclic shift of vertex labels; orbit size divides n.
    """
    orbits = defaultdict(list)
    for s in survivors:
        rep = canonical_orbit_rep(s, n)
        orbits[rep].append(s)
    return orbits


# ============================================================================
# 3.  Output
# ============================================================================

def write_outputs(orbits, n, out_json, out_md):
    """
    Write the orbit-decomposition manifest (JSON) and a human-readable
    summary (Markdown).

    LOGIC:  iterate orbits in canonical-rep lex order; build one record
            per orbit; serialize to disk.

    PHYSICS:  the JSON is consumed by `cascade_kill_orbit_reps_n9.py`
              -- it gets exactly one M per orbit and runs the Laurent
              cascade on it, with the certainty that success lifts to
              the full orbit by cyclic equivariance.
    """
    rep_keys = sorted(orbits.keys())
    rep_records = []
    for orbit_id, rep in enumerate(rep_keys, start=1):
        members = orbits[rep]
        member_idx = 0  # index of one member in the survivors list (any one)
        rec = {
            "orbit_id": orbit_id,
            "size": len(members),
            "representative": [list(c) for c in rep],
            "vertices_used": sorted(vertices_of(rep)),
            "vertices_missed": sorted(missed_vertices(rep, n)),
        }
        rep_records.append(rec)

    with open(out_json, "w") as f:
        json.dump({"n": n, "n_survivors": sum(len(v) for v in orbits.values()),
                   "n_orbits": len(orbits),
                   "orbits": rep_records}, f, indent=2)

    sizes = sorted({len(v) for v in orbits.values()})
    size_dist = defaultdict(int)
    for v in orbits.values():
        size_dist[len(v)] += 1

    lines = []
    lines.append(f"# n={n} step-1-survivor cyclic orbits\n")
    lines.append(f"Group: $\\mathbb{{Z}}_{{{n}}}$ acting by $v \\to v + s "
                 f"\\pmod{{{n}}}$.  Orbit sizes divide {n}.\n")
    n_total = sum(len(v) for v in orbits.values())
    lines.append(f"Total non-tri step-1 survivors: **{n_total}**.")
    lines.append(f"Total orbits: **{len(orbits)}**.")
    lines.append(f"Distinct orbit sizes: {sizes}\n")
    lines.append("Orbit-size distribution:\n")
    lines.append("| Orbit size | # orbits |")
    lines.append("|---|---:|")
    for sz in sorted(size_dist.keys()):
        lines.append(f"| {sz} | {size_dist[sz]} |")
    lines.append("")
    lines.append(
        "## Per-orbit representatives\n"
        "(One canonical rep per orbit; cascade verification on the rep "
        "lifts to the whole orbit by cyclic equivariance.)\n"
    )
    lines.append("| Orbit | Size | Canonical rep | V_missed |")
    lines.append("|---:|---:|---|---|")
    for rec in rep_records:
        rep_tuple = tuple(tuple(c) for c in rec["representative"])
        lines.append(f"| {rec['orbit_id']} | {rec['size']} | "
                     f"`{format_multiset(rep_tuple)}` | "
                     f"{{{','.join(str(v) for v in rec['vertices_missed'])}}} |")
    with open(out_md, "w") as f:
        f.write("\n".join(lines) + "\n")


# ============================================================================
# 4.  Main
# ============================================================================

def main():
    """
    Top-level driver: enumerate, decompose, write outputs.
    """
    n = 9
    print(f"Enumerating step-1 survivors at n={n}...")
    import time
    t0 = time.time()
    survivors = step1_survivors(n)
    print(f"  Found {len(survivors)} non-tri step-1 survivors in "
          f"{time.time()-t0:.1f} s (expected: 1011).")

    print("Computing cyclic-orbit decomposition...")
    orbits = orbit_decomposition(survivors, n)
    sizes = sorted({len(v) for v in orbits.values()})
    print(f"  {len(orbits)} orbits, size distribution: "
          f"{dict((sz, sum(1 for v in orbits.values() if len(v) == sz)) for sz in sizes)}")

    out_json = os.path.join(OUT_DIR, "orbits_n9.json")
    out_md   = os.path.join(OUT_DIR, "orbits_n9.md")
    write_outputs(orbits, n, out_json, out_md)
    print(f"Wrote {out_json} and {out_md}.")


if __name__ == "__main__":
    main()
