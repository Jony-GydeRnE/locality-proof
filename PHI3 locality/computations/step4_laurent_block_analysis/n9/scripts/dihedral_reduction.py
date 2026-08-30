"""
Dihedral (D_9) reduction of the 90 cluster orbits.
================================================================================

WHAT THIS SCRIPT DOES, IN PLAIN LANGUAGE
----------------------------------------
The 90 "cluster orbits" we proved have rank=90 are already grouped under
ROTATIONS of the 9-gon. (That is what "Z_9 orbit" means: rotate the polygon
1/9 turn at a time and keep all the diagrams that look the same.)

But the 9-gon also has REFLECTIONS — flipping it left-to-right. If we mod
out by reflections too (i.e. consider rotations + reflections together,
which is the dihedral group D_9 with 18 elements), some of the 90 orbits
will merge into bigger groups. The result is a SHORTER list — the truly
distinct diagrams when we don't care about flipping.

For each resulting D_9 group we report two numbers:

  ROT COPIES     = how many actual diagrams sit inside the group, under
                   rotation only. (Almost always 9, occasionally less if
                   the diagram is rotation-symmetric.)

  REFL COPIES    = how many EXTRA diagrams enter when we add reflection.
                   * 0 means the diagram already looks the same after
                     flipping (a "palindromic" diagram). The rotation
                     orbit IS the dihedral orbit.
                   * Equal to ROT COPIES means the flipped version is a
                     different rotation orbit, and the two merge into a
                     single D_9 group of size 2 × (rotation orbit size).

INPUTS / OUTPUTS
----------------
INPUTS:
  n9/outputs/orbits_n9.json        -- all 113 cyclic-orbit representatives
  n9/outputs/cluster_partial.json  -- the 90 orbit IDs in the cluster

OUTPUTS:
  n9/outputs/dihedral_reduction.json     -- machine-readable groups + stats
  n9/outputs/dihedral_reduction.md       -- human-readable per-group table
  n9/outputs/cluster_dihedral_gallery.pdf -- one diagram per D_9 group
"""

import json
import math
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Polygon as MplPolygon


# ============================================================================
# 1.  WHERE FILES LIVE
# ============================================================================
N = 9
HERE        = os.path.dirname(os.path.abspath(__file__))
OUTPUTS_DIR = os.path.normpath(os.path.join(HERE, "..", "outputs"))
ORBITS_PATH  = os.path.join(OUTPUTS_DIR, "orbits_n9.json")
CLUSTER_PATH = os.path.join(OUTPUTS_DIR, "cluster_partial.json")
OUT_JSON     = os.path.join(OUTPUTS_DIR, "dihedral_reduction.json")
OUT_MD       = os.path.join(OUTPUTS_DIR, "dihedral_reduction.md")
OUT_PDF      = os.path.join(OUTPUTS_DIR, "cluster_dihedral_gallery.pdf")
SEED_ORBITS  = {22, 46, 88, 108}  # the 4 cascade-failure seeds


# ============================================================================
# 2.  GROUP THEORY -- HOW ROTATIONS AND REFLECTIONS ACT ON 9-GON CHORDS
# ============================================================================
# A 9-gon has vertices labeled 1..9 going clockwise around a circle.
#
# A ROTATION by k notches moves vertex i to vertex i+k (mod 9).
# A REFLECTION across the axis through vertex 1 sends vertex i to (2-i) mod 9.
#
# A chord is an unordered pair {a, b} of non-adjacent vertices. Any
# rotation/reflection acts on a chord by acting on its endpoints and
# re-sorting (so the chord stays "smaller-vertex first"). It acts on a
# multiset of chords by acting on each chord and re-sorting the list.

def rotate_vertex(i, k, n=N):
    """Rotate vertex i by k notches around the n-gon."""
    return ((i - 1 + k) % n) + 1


def reflect_vertex(i, n=N):
    """Reflect vertex i across the axis through vertex 1.
       Sends 1->1, 2<->9, 3<->8, 4<->7, 5<->6."""
    return ((1 - i) % n) + 1


def transform_chord(chord, vertex_perm):
    """Apply vertex permutation to a chord and re-sort to (smaller, larger)."""
    a = vertex_perm(chord[0])
    b = vertex_perm(chord[1])
    return (a, b) if a < b else (b, a)


def transform_multiset(multiset, vertex_perm):
    """Apply vertex permutation to all chords in a multiset, then sort the
       chord list lexicographically -- so the result is in canonical form
       for that particular permutation."""
    return tuple(sorted(transform_chord(c, vertex_perm) for c in multiset))


def canonical_under_rotations(multiset):
    """Return the smallest multiset (lex order) reachable by any rotation.
       Two multisets share this canonical form iff they are in the same
       Z_9 (rotation) orbit."""
    return min(transform_multiset(multiset, lambda i, k=k: rotate_vertex(i, k))
               for k in range(N))


# ============================================================================
# 3.  LOAD MANIFESTS AND BUILD A "MULTISET -> ORBIT ID" LOOKUP
# ============================================================================

def load_inputs():
    """Read the orbit manifest (113 orbits) and the cluster list (90 IDs).
       Return:
         orbit_rep   : dict orbit_id -> chord multiset (as a tuple of tuples)
         orbit_size  : dict orbit_id -> size of the Z_9 orbit (8 or 9 typ.)
         cluster_ids : set of the 90 cluster orbit IDs
    """
    with open(ORBITS_PATH) as f:
        manifest = json.load(f)
    with open(CLUSTER_PATH) as f:
        cluster = json.load(f)

    orbit_rep   = {o["orbit_id"]: tuple(tuple(c) for c in o["representative"])
                   for o in manifest["orbits"]}
    orbit_size  = {o["orbit_id"]: o["size"] for o in manifest["orbits"]}
    cluster_ids = set(cluster["cluster_orbits"])
    return orbit_rep, orbit_size, cluster_ids


# ============================================================================
# 4.  FOR EACH CLUSTER ORBIT, FIND ITS REFLECTION PARTNER ORBIT
# ============================================================================
# The plan: take the orbit's representative, REFLECT it, then ask "which
# Z_9 orbit is this reflected multiset in?"  We answer by reducing the
# reflected multiset to its canonical-under-rotations form and looking it
# up in a precomputed table.

def find_reflection_partners(orbit_rep, cluster_ids):
    """Returns dict mapping each cluster orbit_id to the orbit_id of its
       reflection partner (which may equal itself if the diagram is
       reflection-symmetric, or may be an orbit OUTSIDE the cluster -- we
       record that too)."""
    canonical_to_orbit = {canonical_under_rotations(orbit_rep[oid]): oid
                          for oid in orbit_rep}
    partners = {}
    for oid in cluster_ids:
        M_reflected = transform_multiset(orbit_rep[oid], reflect_vertex)
        canon = canonical_under_rotations(M_reflected)
        partners[oid] = canonical_to_orbit[canon]
    return partners


# ============================================================================
# 5.  GROUP THE CLUSTER ORBITS INTO DIHEDRAL EQUIVALENCE CLASSES
# ============================================================================
# Union-find: two orbits get merged whenever one is the reflection partner
# of the other. We only merge inside the cluster (so a partner sitting in a
# non-cluster orbit gets noted but doesn't add anything to the cluster
# groups -- it's a "leak" out of the cluster under reflection).

def group_into_dihedral_classes(partners, cluster_ids):
    parent = {oid: oid for oid in cluster_ids}

    def find(x):
        # path-compressed union-find
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for oid in cluster_ids:
        partner = partners[oid]
        if partner in cluster_ids:
            union(oid, partner)
        # else: partner is outside the cluster -- record but don't merge.

    groups = {}
    for oid in cluster_ids:
        groups.setdefault(find(oid), []).append(oid)
    # Sort: by smallest orbit ID inside each group.
    return sorted([sorted(grp) for grp in groups.values()],
                  key=lambda g: g[0])


# ============================================================================
# 6.  BUILD A SUMMARY ROW PER D_9 GROUP
# ============================================================================
# For each group we record:
#   - representative orbit ID (smallest in the group)
#   - the chord multiset of that representative
#   - which Z_9 orbits (1 or 2) sit inside the group
#   - rotation copies and reflection copies (see header comment for the
#     plain-language definition)
#   - "kind":
#       palindromic       -- 1 Z_9 orbit, it equals its own reflection orbit
#       palindromic_leak  -- 1 Z_9 orbit, but its reflection lives OUTSIDE
#                            the cluster (so the cluster isn't closed
#                            under D_9 at this orbit)
#       reflection_pair   -- 2 distinct Z_9 orbits, paired by reflection

def summarize_groups(groups, partners, orbit_rep, orbit_size, cluster_ids):
    summary = []
    for grp in groups:
        rep_oid = grp[0]
        rot_copies = sum(orbit_size[o] for o in grp)
        if len(grp) == 1:
            # The orbit's reflection might fall inside or outside the cluster.
            partner_outside = (partners[rep_oid] not in cluster_ids)
            if partner_outside:
                kind = "palindromic_leak"
                refl_copies = 0
                external_partner = partners[rep_oid]
            else:
                kind = "palindromic"
                refl_copies = 0
                external_partner = None
        elif len(grp) == 2:
            kind = "reflection_pair"
            refl_copies = rot_copies   # equal-size second half
            external_partner = None
        else:
            # Mathematically a D_9 group of cluster orbits has at most 2 Z_9
            # orbits inside it. If we ever see > 2, something is wrong.
            kind = "unexpected_larger_group"
            refl_copies = -1
            external_partner = None
        summary.append({
            "rep_orbit_id": rep_oid,
            "rep_multiset": [list(c) for c in orbit_rep[rep_oid]],
            "z9_orbits_in_group": grp,
            "n_z9_orbits": len(grp),
            "rot_copies": rot_copies,
            "refl_copies": refl_copies,
            "kind": kind,
            "external_partner_orbit": external_partner,
            "is_seed": rep_oid in SEED_ORBITS or any(o in SEED_ORBITS for o in grp),
        })
    return summary


# ============================================================================
# 7.  DRAW DIAGRAMS (same style as the existing 90-orbit gallery)
# ============================================================================

def vertex_xy(i, n=N, R=1.0):
    """Cartesian coords of vertex i on a regular n-gon (vertex 1 at top)."""
    angle = math.pi / 2 - 2 * math.pi * (i - 1) / n
    return (R * math.cos(angle), R * math.sin(angle))


def draw_diagram(ax, multiset, n=N):
    """Render: filled black 9-gon with red chord overlays for the multiset.
       Multiplicity-2 chords get drawn as two parallel red lines."""
    coords = [vertex_xy(i, n) for i in range(1, n + 1)]
    ax.add_patch(MplPolygon(coords, closed=True,
                            facecolor='black', edgecolor='black'))
    counts = {}
    for c in multiset:
        counts[tuple(c)] = counts.get(tuple(c), 0) + 1
    for chord, mult in counts.items():
        i, j = chord
        x1, y1 = vertex_xy(i)
        x2, y2 = vertex_xy(j)
        for k in range(mult):
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy) or 1.0
            nx, ny = -dy / length, dx / length
            shift = 0.03 * (k - (mult - 1) / 2.0)
            ax.plot([x1 + shift * nx, x2 + shift * nx],
                    [y1 + shift * ny, y2 + shift * ny],
                    color='red', linewidth=2.4)
    for v in range(1, n + 1):
        x, y = vertex_xy(v, R=1.12)
        ax.text(x, y, str(v), fontsize=6, ha='center', va='center')
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect('equal')
    ax.axis('off')


def caption_for_group(grp_summary):
    """One-paragraph caption for a D_9 group's diagram."""
    z9 = grp_summary["z9_orbits_in_group"]
    seed_marker = "  [SEED]" if grp_summary["is_seed"] else ""
    if grp_summary["kind"] == "palindromic":
        kind_str = f"palindromic (rot copies={grp_summary['rot_copies']}, no reflection siblings)"
    elif grp_summary["kind"] == "palindromic_leak":
        kind_str = (f"palindromic-leak (Z_9-only inside cluster; reflection "
                    f"partner = orbit {grp_summary['external_partner_orbit']} "
                    f"is OUTSIDE the cluster)")
    elif grp_summary["kind"] == "reflection_pair":
        kind_str = (f"reflection pair: Z_9 orbits {z9[0]} and {z9[1]} "
                    f"(rot copies={grp_summary['rot_copies']}, refl copies={grp_summary['refl_copies']})")
    else:
        kind_str = grp_summary["kind"]
    M_str = '{' + ', '.join(f"({c[0]},{c[1]})" for c in grp_summary["rep_multiset"]) + '}'
    return f"D_9 group rep = orbit #{grp_summary['rep_orbit_id']}{seed_marker}\nZ_9 ⊂ : {z9}   {kind_str}\n{M_str}"


def render_pdf(summary):
    """Render one page per ~6 D_9 groups, with a cover page summary."""
    n_groups = len(summary)
    per_page = 6
    cols = 2
    rows = per_page // cols
    n_pages = math.ceil(n_groups / per_page)

    n_palindrome = sum(1 for g in summary if g["kind"] == "palindromic")
    n_pleak      = sum(1 for g in summary if g["kind"] == "palindromic_leak")
    n_pair       = sum(1 for g in summary if g["kind"] == "reflection_pair")

    with PdfPages(OUT_PDF) as pdf:
        # ----- cover page -----
        cover = plt.figure(figsize=(8.5, 11))
        cover.text(0.5, 0.85, "n=9 cluster orbits, modded out by reflections",
                   ha='center', fontsize=18, weight='bold')
        cover.text(0.5, 0.79, f"D_9 (rotations + reflections) reduction",
                   ha='center', fontsize=12)
        cover.text(0.5, 0.70,
                   f"Original cluster: 90 cyclic orbits (Z_9)\n"
                   f"After modding out by reflection: {n_groups} D_9 groups\n",
                   ha='center', fontsize=11)
        cover.text(0.5, 0.58,
                   f"Breakdown of the {n_groups} groups:\n"
                   f"  - palindromic (reflection-self-symmetric, alone in group): "
                   f"{n_palindrome}\n"
                   f"  - reflection pairs (two Z_9 orbits merged into one): "
                   f"{n_pair}\n"
                   f"  - palindromic with reflection partner OUTSIDE cluster: "
                   f"{n_pleak}",
                   ha='center', fontsize=10)
        cover.text(0.5, 0.42,
                   "Read each diagram caption like this:\n"
                   "  ‘rep = orbit #X’  →  the smallest Z_9 orbit ID in the group\n"
                   "  ‘Z_9 ⊂ : [a, b]’  →  which Z_9 orbits sit inside the D_9 group\n"
                   "  ‘rot copies = N’  →  total number of distinct diagrams under rotation only\n"
                   "  ‘refl copies = N’ →  EXTRA diagrams that show up only when reflection is added",
                   ha='center', fontsize=9, style='italic')
        cover.text(0.5, 0.06,
                   f"Pages: {n_pages}, 6 diagrams/page",
                   ha='center', fontsize=8)
        pdf.savefig(cover); plt.close(cover)

        # ----- diagram pages -----
        for page_idx in range(n_pages):
            start = page_idx * per_page
            page = summary[start: start + per_page]
            fig, axes = plt.subplots(rows, cols, figsize=(8.5, 11))
            fig.suptitle(
                f"n=9 D_9 cluster groups  —  page {page_idx + 1} / {n_pages}",
                fontsize=11)
            axes_flat = axes.flatten()
            for ax in axes_flat:
                ax.axis('off')
            for ax, grp in zip(axes_flat, page):
                draw_diagram(ax, grp["rep_multiset"])
                ax.set_title(caption_for_group(grp), fontsize=7, pad=4)
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            pdf.savefig(fig); plt.close(fig)
            print(f"  page {page_idx + 1}/{n_pages} done")
    print(f"Wrote: {OUT_PDF}")


# ============================================================================
# 8.  WRITE THE JSON + MARKDOWN SUMMARY
# ============================================================================

def write_text_outputs(summary, partners, cluster_ids):
    n_groups = len(summary)
    n_palindrome = sum(1 for g in summary if g["kind"] == "palindromic")
    n_pleak      = sum(1 for g in summary if g["kind"] == "palindromic_leak")
    n_pair       = sum(1 for g in summary if g["kind"] == "reflection_pair")

    with open(OUT_JSON, "w") as f:
        json.dump({
            "n_z9_cluster_orbits": len(cluster_ids),
            "n_d9_groups": n_groups,
            "n_palindromic": n_palindrome,
            "n_palindromic_with_external_partner": n_pleak,
            "n_reflection_pairs": n_pair,
            "groups": summary,
            "reflection_partner_of_each_cluster_orbit": partners,
        }, f, indent=2)

    lines = []
    lines.append("# D_9 reduction of the 90 n=9 cluster orbits\n")
    lines.append(f"- Original cluster (Z_9 orbits): **{len(cluster_ids)}**")
    lines.append(f"- After modding out by reflection (D_9): **{n_groups}** groups\n")
    lines.append("## Breakdown\n")
    lines.append("| Kind | Count | Meaning |")
    lines.append("|---|---:|---|")
    lines.append(f"| Palindromic | {n_palindrome} | One Z_9 orbit, equals its own reflection orbit. |")
    lines.append(f"| Reflection pair | {n_pair} | Two Z_9 orbits, swapped by reflection. |")
    lines.append(f"| Palindromic-leak | {n_pleak} | One Z_9 orbit inside cluster; reflection partner is OUTSIDE the cluster. |")
    lines.append("")
    lines.append(f"Total Z_9 orbits accounted for: "
                 f"{n_palindrome + n_pleak + 2*n_pair} (should equal "
                 f"{len(cluster_ids)}).\n")
    lines.append("## Per-group table\n")
    lines.append("| D_9 group # | Rep orbit | Z_9 orbits in group | Kind | Rot copies | Refl copies | Multiset |")
    lines.append("|---:|---:|---|---|---:|---:|---|")
    for i, g in enumerate(summary, 1):
        ms = '{' + ', '.join(f"({c[0]},{c[1]})" for c in g["rep_multiset"]) + '}'
        seed_tag = " ★" if g["is_seed"] else ""
        z9_str = ", ".join(str(o) for o in g["z9_orbits_in_group"])
        lines.append(f"| {i}{seed_tag} | {g['rep_orbit_id']} | {z9_str} | "
                     f"{g['kind']} | {g['rot_copies']} | {g['refl_copies']} | "
                     f"`{ms}` |")
    lines.append("")
    lines.append("★ = the D_9 group contains one of the 4 cascade-failure seed "
                 "orbits {22, 46, 88, 108}.")
    with open(OUT_MD, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote: {OUT_JSON}")
    print(f"Wrote: {OUT_MD}")


# ============================================================================
# 9.  MAIN
# ============================================================================

def main():
    orbit_rep, orbit_size, cluster_ids = load_inputs()
    print(f"Cluster: {len(cluster_ids)} Z_9 orbits.")

    partners = find_reflection_partners(orbit_rep, cluster_ids)
    groups = group_into_dihedral_classes(partners, cluster_ids)
    summary = summarize_groups(groups, partners, orbit_rep,
                               orbit_size, cluster_ids)

    print(f"After D_9 reduction: {len(groups)} groups.")
    print("Breakdown:")
    print(f"  palindromic                          : "
          f"{sum(1 for g in summary if g['kind']=='palindromic')}")
    print(f"  reflection pair                      : "
          f"{sum(1 for g in summary if g['kind']=='reflection_pair')}")
    print(f"  palindromic with external partner    : "
          f"{sum(1 for g in summary if g['kind']=='palindromic_leak')}")

    write_text_outputs(summary, partners, cluster_ids)
    render_pdf(summary)

    try:
        if sys.platform == "darwin":
            subprocess.Popen(["open", OUT_PDF])
        elif sys.platform.startswith("linux"):
            subprocess.Popen(["xdg-open", OUT_PDF])
        elif sys.platform == "win32":
            os.startfile(OUT_PDF)  # type: ignore[attr-defined]
        print("Opened gallery PDF.")
    except Exception as e:
        print(f"(Could not auto-open: {e})")


if __name__ == "__main__":
    main()
