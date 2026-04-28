"""
================================================================================
visualize_untouched.py  --  Print the 4 untouched Step-2 components at n=9
                              explicitly, with diagrams and Step-2 equations
================================================================================

PURPOSE
-------
After bridging TEST 1 (depth-1 from 23 non-cluster orbits), the
n=9 Step-2 flip-graph still has 4 components that NO step-1 / step-3
mechanism touches: components {1, 3, 7, 14}.  This script prints them
out for visual inspection, including:

  - the orbit IDs in each component,
  - one canonical-representative triangulation per orbit (as a set of
    6 chords on the 9-gon),
  - an ASCII diagram of the 9-gon with chords drawn for one
    representative per component,
  - the Step-2 swap equations *within* the component (orbit-level
    edges that link the orbit's triangulations among themselves),
  - structural features (cyclic stabilizer order, chord-length
    distribution, ear vertices),
  - if matplotlib is available, also save PNG diagrams to diagrams/.

Output:  outputs/untouched_components_readable.md
         (and diagrams/component_<id>_orbit_<id>.png if matplotlib OK)

CODE STYLE: per repo standards, every function has LOGIC + PHYSICS
docstring.
"""

import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
N8_DIR = os.path.normpath(os.path.join(
    HERE, "..", "..", "step4_laurent_block_analysis", "n8", "scripts"))
sys.path.insert(0, N8_DIR)
from cascade_kill_n8 import format_multiset  # noqa: E402
from analyze_recipes import shift_multiset  # noqa: E402

OUT_MD = os.path.join(HERE, "outputs", "untouched_components_readable.md")
DIAG_DIR = os.path.join(HERE, "diagrams")


def load_components_data():
    """Load orbit reps + Step-2 components + edges."""
    with open(os.path.join(HERE, "outputs",
                           "n9_triangulation_orbits.json")) as f:
        tri_data = json.load(f)
    rep_to_tri_id = {}
    tri_id_to_rep = {}
    for o in tri_data["orbits"]:
        rep = tuple(tuple(c) for c in o["representative"])
        rep_to_tri_id[rep] = o["orbit_id"]
        tri_id_to_rep[o["orbit_id"]] = rep
        tri_id_to_rep_full = tri_id_to_rep
    with open(os.path.join(HERE, "outputs",
                           "n9_step2_bare_swap_pairs.json")) as f:
        swap_data = json.load(f)
    edges = [(e["orbit_pair"][0], e["orbit_pair"][1])
             for e in swap_data["orbit_edges"]]
    # Edges with full info
    edges_full = swap_data["orbit_edges"]
    # Triangulation-level edges
    tri_edges = swap_data["triangulation_edges"]
    return tri_id_to_rep, edges, edges_full, tri_edges


def union_find_components(nodes, edges):
    """Compute connected components."""
    parent = {n: n for n in nodes}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    comps = defaultdict(list)
    for n in nodes:
        comps[find(n)].append(n)
    # sort components by min orbit id, then assign component IDs in that order
    sorted_roots = sorted(comps.keys(), key=lambda r: min(comps[r]))
    cid_map = {r: i + 1 for i, r in enumerate(sorted_roots)}
    return {n: cid_map[find(n)] for n in nodes}, cid_map


def cyclic_stabilizer(rep, n=9):
    """Return the order of the cyclic stabilizer subgroup of rep."""
    for s in range(1, n + 1):
        if shift_multiset(rep, s, n) == tuple(sorted(rep)):
            return n // s
    return 1


def chord_length(c, n=9):
    """Cyclic distance between the two endpoints of chord c."""
    a, b = c
    d = abs(b - a)
    return min(d, n - d)


def ear_vertices(rep, n=9):
    used = set()
    for (i, j) in rep:
        used.add(i); used.add(j)
    return sorted(set(range(1, n + 1)) - used)


def ascii_polygon_diagram(rep, n=9, label=""):
    """
    Build a small ASCII polygon diagram with chords listed.

    LOGIC:  produce a 9-vertex circular layout indicator + the chord
            list with cyclic-distance labels.

    PHYSICS:  small visual aid; the 9-gon is too dense for a useful
              ASCII picture, so we list the chords with their lengths
              and ears.  PNG diagrams via matplotlib are the better
              visualisation.
    """
    lines = []
    lines.append(f"### {label}")
    lines.append("```")
    lines.append("9-gon vertex layout (cyclic):")
    lines.append("        1")
    lines.append("    9       2")
    lines.append("  8           3")
    lines.append("  7           4")
    lines.append("    6       5")
    lines.append("        ")
    lines.append("Chords in T (with cyclic length):")
    for (i, j) in sorted(rep):
        L = chord_length((i, j))
        lines.append(f"    ({i},{j})  length = {L}")
    lines.append("Ear vertices (no chord touches): " + str(ear_vertices(rep)))
    lines.append("```")
    return "\n".join(lines)


def make_png_diagram(rep, n=9, save_path=None):
    """
    Try to draw the 9-gon with chords using matplotlib.

    LOGIC:  vertices on a regular 9-gon at unit circle; chords as line
            segments.  Save PNG if matplotlib is available.

    PHYSICS:  cleaner visualisation than ASCII.  The user can scan
              chord patterns visually to spot common features across
              the 4 untouched components.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return False
    fig, ax = plt.subplots(figsize=(6, 6))
    angles = np.array([2 * np.pi * (k - 1) / n + np.pi / 2
                        for k in range(1, n + 1)])
    xs = np.cos(angles)
    ys = np.sin(angles)
    # Polygon edges
    for i in range(n):
        nxt = (i + 1) % n
        ax.plot([xs[i], xs[nxt]], [ys[i], ys[nxt]],
                color='lightgray', linewidth=1)
    # Vertices + labels
    ax.scatter(xs, ys, color='black', zorder=3)
    for k in range(1, n + 1):
        ax.text(xs[k-1] * 1.10, ys[k-1] * 1.10, str(k),
                ha='center', va='center', fontsize=12)
    # Chords
    for (i, j) in rep:
        ax.plot([xs[i-1], xs[j-1]], [ys[i-1], ys[j-1]],
                color='blue', linewidth=2)
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect('equal')
    ax.axis('off')
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=120, bbox_inches='tight')
    plt.close(fig)
    return True


def main():
    os.makedirs(os.path.dirname(OUT_MD), exist_ok=True)
    os.makedirs(DIAG_DIR, exist_ok=True)

    tri_id_to_rep, edges, edges_full, tri_edges = load_components_data()
    nodes = list(tri_id_to_rep.keys())
    node_to_cid, cid_map = union_find_components(nodes, edges)
    # Group orbits by component
    components = defaultdict(list)
    for oid, cid in node_to_cid.items():
        components[cid].append(oid)

    untouched_target = {1, 3, 7, 14}

    lines = []
    lines.append("# n=9 Step-2 untouched triangulation components — readable\n")
    lines.append("This file prints the four Step-2 components that the "
                 "depth-1 fingerprint matrix on the 90-orbit cluster does "
                 "NOT bridge, and that the depth-1 fingerprints of the 23 "
                 "non-cluster non-tri orbits do not bridge either "
                 "(TEST 1 result).\n")
    lines.append(f"Components: **{sorted(untouched_target)}**\n")
    lines.append(f"Each component below lists its triangulation orbits "
                 f"(canonical reps), Step-2 swap edges within the "
                 f"component, structural features, and ASCII / PNG "
                 f"diagrams of one representative per orbit.\n")

    # Build orbit-level edges with full info for printing
    orbit_edges_full = defaultdict(list)
    for e in edges_full:
        a, b = e["orbit_pair"]
        orbit_edges_full[(a, b)].append(e)

    # Build triangulation-level edges for printing within-component swap eqs
    tri_edges_by_orbit_pair = defaultdict(list)
    for te in tri_edges:
        key = tuple(sorted([te["T1_orbit"], te["T2_orbit"]]))
        tri_edges_by_orbit_pair[key].append(te)

    for cid in sorted(untouched_target):
        orbit_ids = sorted(components[cid])
        lines.append(f"\n## Component {cid}\n")
        lines.append(f"Orbits in component: **{orbit_ids}** "
                     f"({len(orbit_ids)} orbits)\n")

        # Per-orbit representatives + structural features
        lines.append("### Triangulation orbit representatives\n")
        lines.append("| Orbit ID | Rep (chord set) | Stabilizer order | "
                     "Chord lengths | Ear vertices |")
        lines.append("|---:|---|---:|---|---|")
        for oid in orbit_ids:
            rep = tri_id_to_rep[oid]
            stab = cyclic_stabilizer(rep, 9)
            lens = sorted(chord_length(c) for c in rep)
            ears = ear_vertices(rep, 9)
            chord_str = ", ".join(f"({i},{j})" for (i, j) in sorted(rep))
            lines.append(f"| {oid} | {{{chord_str}}} | {stab} | "
                         f"{lens} | {ears} |")
        lines.append("")

        # Step-2 swap edges within the component
        lines.append(f"### Step-2 swap edges *within* component {cid}\n")
        in_comp_edges = []
        for (a, b), edges_for_pair in orbit_edges_full.items():
            if a in orbit_ids and b in orbit_ids:
                in_comp_edges.append((a, b, edges_for_pair[0]))
        if not in_comp_edges:
            lines.append("(none)")
        else:
            lines.append("Each row gives one orbit-level Step-2 edge: "
                         "two triangulation orbits whose representatives "
                         "are related by a single bare/special swap at "
                         "some cyclic zone(s). Any two orbits in this "
                         "component are linked by a chain of such edges.\n")
            lines.append("| Orbit A | Orbit B | Sample triangulation pair "
                         "(at zones) |")
            lines.append("|---:|---:|---|")
            for a, b, e in in_comp_edges:
                # Find one triangulation pair from this orbit pair
                key = tuple(sorted([a, b]))
                te_list = tri_edges_by_orbit_pair.get(key, [])
                if te_list:
                    te = te_list[0]
                    T1 = [tuple(c) for c in te["T1"]]
                    T2 = [tuple(c) for c in te["T2"]]
                    zones = te["zones"]
                    lines.append(f"| {a} | {b} | "
                                 f"`{format_multiset(tuple(T1))}` ↔ "
                                 f"`{format_multiset(tuple(T2))}` "
                                 f"(zones {zones}) |")
                else:
                    lines.append(f"| {a} | {b} | (no triangulation-level "
                                 f"detail recorded) |")
        lines.append("")

        # Diagrams + ASCII per orbit
        lines.append("### Diagrams\n")
        for oid in orbit_ids:
            rep = tri_id_to_rep[oid]
            lines.append("")
            lines.append(ascii_polygon_diagram(
                rep, 9, label=f"Orbit {oid}"))
            png_path = os.path.join(
                DIAG_DIR, f"component_{cid}_orbit_{oid}.png")
            ok = make_png_diagram(rep, 9, save_path=png_path)
            rel = os.path.relpath(png_path, os.path.dirname(OUT_MD))
            if ok:
                lines.append(f"\n![Component {cid} orbit {oid}]({rel})\n")
        lines.append("")

    # Final summary section
    lines.append("\n---\n")
    lines.append("## Why these components are isolated\n")
    lines.append(
        "**Step 1 (layer-0 kill)** never touches triangulations at all "
        "(local-survival lemma: no triangulation is killable at any zone). "
        "So Step-1 contributes zero equations involving these orbits.\n")
    lines.append(
        "**Step 2 (bare/special swap on triangulations)** equates "
        "coefficients within each Step-2 component. By construction the "
        "4 components above are CLOSED under this swap — every "
        "bare/special swap from any of their triangulations either "
        "stays in the same component or fails (the swapped multiset "
        "isn't a triangulation).\n")
    lines.append(
        "**Step 3 (block-rule cluster matrix at depth-1)** equates "
        "cluster non-tri coefficients to triangulation cousins via "
        "fingerprint equations. The 32 triangulation cousins from the "
        "cluster matrix land in 12 of the 16 components; these 4 "
        "components have NO triangulation cousin coupled to any of "
        "the 90 cluster orbits at depth-1.\n")
    lines.append(
        "**TEST 1 (depth-1 fingerprints from 23 non-cluster non-tri "
        "orbits)** also yields zero new bridges into these components.\n")
    lines.append(
        "Combined: every depth-1 hidden-zero relation either has all "
        "its triangulation cousins in the *other* 12 components, or "
        "has no triangulation cousins at all. The 4 untouched "
        "components are *invisible* to the depth-1 mechanism.\n")
    lines.append("## What might bridge them\n")
    lines.append(
        "Possibilities for the missing mechanism, all using ONLY "
        "1-zero hidden-zero conditions:\n"
        "1. **Depth-2 (or deeper) Laurent fingerprints** that pull in "
        "triangulation cousins these depth-1 equations miss. (Test 2 "
        "in `build_bridges.py`.)\n"
        "2. **Pure-triangulation Laurent extractions** — free-chord "
        "monomials at low order whose only contributors are "
        "triangulations, giving triangulation-only relations.\n"
        "3. **Multi-bare Step-2 identities** — currently the swap "
        "exchanges one bare for one special; allowing multiple "
        "bare/special exchanges in a single multiset may produce new "
        "edges.\n")

    with open(OUT_MD, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Wrote {OUT_MD}")
    print(f"Diagrams (if matplotlib available): {DIAG_DIR}/")


if __name__ == "__main__":
    main()
