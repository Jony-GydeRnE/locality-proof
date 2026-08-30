"""
Length-spectrum reduction: collapse the 55 D_9 cluster groups by chord-length type.
================================================================================

WHAT THIS SCRIPT IS FOR
-----------------------
The previous reduction step took the 90 cyclic-orbit (Z_9) cluster
representatives and grouped them under reflections, leaving 55 D_9
groups. But many of those 55 still have visually-similar "shapes" —
they share the same SPECTRUM OF CHORD LENGTHS, even when they aren't
related by any rotation or reflection of the 9-gon.

This script applies a coarser equivalence:

       ┌─────────────────────────────────────────────────────────────┐
       │  Two multisets are SHAPE-EQUIVALENT iff their sorted         │
       │  chord-length multisets are equal.                           │
       └─────────────────────────────────────────────────────────────┘

A chord {a, b} in the 9-gon has LENGTH = the shorter cyclic distance
between vertices a and b, so length ∈ {2, 3, 4}. (Length 1 would be a
polygon edge, which we don't allow as a "chord". Length 4 is the
"longest" chord on a 9-gon.) Each multiset has 6 chords, so its
length-spectrum is a 6-element multiset of values from {2, 3, 4}.

WHY THIS IS A COARSER EQUIVALENCE THAN D_9
------------------------------------------
The dihedral group D_9 acts by rotation and reflection of the 9-gon.
Both preserve the cyclic distance between any two vertices, so D_9
preserves chord lengths. Therefore:

     two multisets D_9-equivalent  ⇒  same length-spectrum

But the converse FAILS: there exist multisets with the same
length-spectrum that are not D_9-equivalent (they live on different
"slots" around the polygon). So:

  same length-spectrum  ⊋  same D_9 orbit

i.e. the length-spectrum equivalence is strictly coarser. It groups
together every multiset that "uses the same kinds of chords", even when
those chords live in different geometric positions.

WHAT YOU SHOULD READ THIS AS
----------------------------
The 90 Z_9 cluster orbits collapse to 55 D_9 groups (rotations +
reflections) and then collapse FURTHER to a much smaller number of
length-spectrum shapes. The shape number is the count of "kinds of
6-chord configurations" the cluster sees, ignoring all positional
information.

INPUTS:
  n9/outputs/dihedral_reduction.json   -- 55 D_9 groups + reps
  n9/outputs/orbits_n9.json            -- chord multisets

OUTPUTS:
  n9/outputs/length_spectrum_reduction.json
  n9/outputs/length_spectrum_reduction.md
  n9/outputs/length_spectrum_gallery.pdf  -- one diagram per shape +
                                             cover page explaining the rule
"""

import json
import math
import os
import subprocess
import sys
from collections import defaultdict, Counter

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
DIHEDRAL_PATH = os.path.join(OUTPUTS_DIR, "dihedral_reduction.json")
ORBITS_PATH   = os.path.join(OUTPUTS_DIR, "orbits_n9.json")
OUT_JSON = os.path.join(OUTPUTS_DIR, "length_spectrum_reduction.json")
OUT_MD   = os.path.join(OUTPUTS_DIR, "length_spectrum_reduction.md")
OUT_PDF  = os.path.join(OUTPUTS_DIR, "length_spectrum_gallery.pdf")
SEED_ORBITS = {22, 46, 88, 108}


# ============================================================================
# 2.  CHORD LENGTH AND LENGTH-SPECTRUM
# ============================================================================
# A chord {a, b} on the 9-gon: its LENGTH is the shorter cyclic distance
# between vertices a and b, so length = min(|a-b|, 9-|a-b|), which is 2,
# 3, or 4. (Length 1 = polygon edge, not a chord.)

def chord_length(chord, n=N):
    a, b = chord
    d = abs(b - a)
    return min(d, n - d)


def length_spectrum(multiset, n=N):
    """Return the sorted length-multiset of a chord multiset, as a tuple
       of integers (so it's hashable for grouping)."""
    return tuple(sorted(chord_length(c, n) for c in multiset))


def spectrum_summary(spec):
    """Pretty-print like '2^5 · 3^1' so it's easy to read."""
    cnt = Counter(spec)
    parts = []
    for L in sorted(cnt):
        parts.append(f"{L}^{cnt[L]}" if cnt[L] > 1 else str(L))
    return " · ".join(parts)


# ============================================================================
# 3.  LOAD INPUTS AND GROUP BY SPECTRUM
# ============================================================================

def load_inputs():
    with open(DIHEDRAL_PATH) as f:
        d9 = json.load(f)
    with open(ORBITS_PATH) as f:
        orbs = json.load(f)
    orbit_size = {o["orbit_id"]: o["size"] for o in orbs["orbits"]}
    return d9, orbit_size


def group_by_spectrum(d9_data, orbit_size):
    """Return dict: spectrum_tuple -> list of D_9 group dicts that have
       that spectrum, plus a 'rep_multiset' on each.
       Also augment each D_9 group with computed totals."""
    groups_by_spec = defaultdict(list)
    for grp in d9_data["groups"]:
        ms = [tuple(c) for c in grp["rep_multiset"]]
        spec = length_spectrum(ms)
        # Compute total Z_9-orbit-element count contained in this D_9 group
        # (it's the sum of orbit sizes of all Z_9 orbits in the group).
        total_multisets = sum(orbit_size[o] for o in grp["z9_orbits_in_group"])
        groups_by_spec[spec].append({
            **grp,
            "rep_multiset_tuples": ms,
            "total_multisets_in_d9_group": total_multisets,
        })
    return groups_by_spec


# ============================================================================
# 4.  DRAWING (same primitives as the other galleries)
# ============================================================================

def vertex_xy(i, n=N, R=1.0):
    angle = math.pi / 2 - 2 * math.pi * (i - 1) / n
    return (R * math.cos(angle), R * math.sin(angle))


def draw_diagram(ax, multiset, n=N, highlight_doubles=True):
    coords = [vertex_xy(i, n) for i in range(1, n + 1)]
    ax.add_patch(MplPolygon(coords, closed=True,
                            facecolor='black', edgecolor='black'))
    counts = Counter(tuple(c) for c in multiset)
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


# ============================================================================
# 5.  RENDER THE FULL PDF
# ============================================================================

def render_pdf(groups_by_spec, n_d9, n_z9, total_multisets):
    n_shapes = len(groups_by_spec)

    # Sort shapes: first by spectrum (lexicographically), then within shape
    # by smallest D_9 group ID.
    sorted_specs = sorted(groups_by_spec.keys())

    # Stats per shape
    shape_stats = []
    for spec in sorted_specs:
        grps = sorted(groups_by_spec[spec],
                      key=lambda g: g["rep_orbit_id"])
        total_z9 = sum(len(g["z9_orbits_in_group"]) for g in grps)
        total_ms = sum(g["total_multisets_in_d9_group"] for g in grps)
        shape_stats.append({
            "spectrum": spec,
            "n_d9_groups": len(grps),
            "n_z9_orbits": total_z9,
            "n_multisets": total_ms,
            "groups": grps,
        })

    with PdfPages(OUT_PDF) as pdf:
        # --------------------------------------------------------------
        # PAGE 1 -- COVER + EXPLANATION OF THE EQUIVALENCE
        # --------------------------------------------------------------
        cov = plt.figure(figsize=(8.5, 11))

        cov.text(0.5, 0.94, "n=9 cluster: chord-length-spectrum reduction",
                 ha='center', fontsize=18, weight='bold')

        cov.text(0.5, 0.89,
                 "(coarser than D_9: groups multisets that 'use the same kinds of chords')",
                 ha='center', fontsize=10, style='italic')

        # The reduction chain
        cov.text(0.5, 0.81, "REDUCTION CHAIN", ha='center', fontsize=12,
                 weight='bold')
        cov.text(0.5, 0.74,
                 f"  1 011  step-1 non-tri survivors\n"
                 f"   ↓  group by cyclic rotation Z_9\n"
                 f"   113 Z_9 orbits  (90 in cluster)\n"
                 f"   ↓  add reflections → D_9\n"
                 f"   55 D_9 groups\n"
                 f"   ↓  forget positions, keep only chord-length spectrum\n"
                 f"   {n_shapes} length-spectrum shapes  ← THIS PDF",
                 ha='center', fontsize=11, family='monospace')

        # The rule
        cov.text(0.5, 0.55, "THE EQUIVALENCE RULE",
                 ha='center', fontsize=12, weight='bold')
        cov.text(0.5, 0.43,
                 "A chord {a, b} on the 9-gon has LENGTH = min(|a-b|, 9-|a-b|),\n"
                 "so length ∈ {2, 3, 4}.  (Length 1 = polygon edge, not a chord.)\n"
                 "\n"
                 "The LENGTH SPECTRUM of a 6-chord multiset is the sorted multiset\n"
                 "of its 6 chord lengths — e.g. (2,2,2,2,2,3) means 'five length-2\n"
                 "chords plus one length-3 chord'.\n"
                 "\n"
                 "Two multisets are SHAPE-EQUIVALENT iff their length spectra match.\n"
                 "\n"
                 "This is COARSER than D_9: rotations and reflections preserve chord\n"
                 "length, but two configurations can have identical length spectra\n"
                 "while sitting in different D_9 orbits (different positional 'slots'\n"
                 "around the polygon).",
                 ha='center', fontsize=10)

        # Headline numbers
        cov.text(0.5, 0.18, "HEADLINE NUMBERS",
                 ha='center', fontsize=12, weight='bold')
        cov.text(0.5, 0.10,
                 f"  Z_9 orbits in cluster   : {n_z9}\n"
                 f"  D_9 groups in cluster   : {n_d9}\n"
                 f"  Length-spectrum shapes  : {n_shapes}    ← total kinds of 'chord-recipe'\n"
                 f"  Total multisets covered : {total_multisets}",
                 ha='center', fontsize=11, family='monospace')

        pdf.savefig(cov); plt.close(cov)

        # --------------------------------------------------------------
        # PAGE 2 -- STATS TABLE BY SHAPE
        # --------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(8.5, 11))
        ax.axis('off')
        ax.set_title(f"All {n_shapes} length-spectrum shapes (sorted by spectrum)",
                     fontsize=13, weight='bold', pad=20)

        cell_text = []
        for i, ss in enumerate(shape_stats, 1):
            spec_str = spectrum_summary(ss["spectrum"])
            cell_text.append([
                str(i),
                spec_str,
                str(ss["n_d9_groups"]),
                str(ss["n_z9_orbits"]),
                str(ss["n_multisets"]),
            ])
        table = ax.table(
            cellText=cell_text,
            colLabels=["shape #", "length spectrum", "# D_9 groups",
                       "# Z_9 orbits", "# multisets"],
            loc='upper center', cellLoc='center',
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.3)
        pdf.savefig(fig); plt.close(fig)

        # --------------------------------------------------------------
        # PAGES 3+ : ONE PAGE PER SHAPE (with diagram + member list)
        # --------------------------------------------------------------
        for i, ss in enumerate(shape_stats, 1):
            spec = ss["spectrum"]
            grps = ss["groups"]
            spec_str = spectrum_summary(spec)
            seed_d9 = [g for g in grps
                       if any(o in SEED_ORBITS for o in g["z9_orbits_in_group"])]

            # Pick the "rep" diagram: smallest D_9 group's rep
            rep_grp = grps[0]
            rep_ms  = rep_grp["rep_multiset_tuples"]

            fig = plt.figure(figsize=(8.5, 11))
            fig.suptitle(
                f"Shape #{i}/{n_shapes}    spectrum  {spec_str}\n"
                f"{ss['n_d9_groups']} D_9 group(s)  ·  "
                f"{ss['n_z9_orbits']} Z_9 orbit(s)  ·  "
                f"{ss['n_multisets']} multisets",
                fontsize=12)

            # representative diagram (top half)
            ax1 = fig.add_axes([0.27, 0.50, 0.46, 0.40])
            draw_diagram(ax1, rep_ms)
            rep_str = '{' + ', '.join(f"({c[0]},{c[1]})" for c in rep_ms) + '}'
            ax1.set_title(f"representative: D_9 grp orbit #{rep_grp['rep_orbit_id']}\n"
                          f"{rep_str}",
                          fontsize=9)

            # text below: list of all D_9 groups with this shape
            text_lines = ["Members of this shape (D_9 group → Z_9 orbits):", ""]
            for g in grps:
                seed_tag = " ★" if any(o in SEED_ORBITS
                                         for o in g["z9_orbits_in_group"]) else ""
                z9_str = ", ".join(str(o) for o in g["z9_orbits_in_group"])
                ms = g["rep_multiset_tuples"]
                ms_short = '{' + ', '.join(f"({c[0]},{c[1]})" for c in ms) + '}'
                text_lines.append(
                    f"  D_9 #{g['rep_orbit_id']:>3}{seed_tag}  "
                    f"({g['kind']:<24}) "
                    f"Z_9: [{z9_str}]   {ms_short}"
                )
            if seed_d9:
                text_lines.append("")
                text_lines.append("  ★ = D_9 group contains one of the cascade-failure "
                                  "seed orbits {22, 46, 88, 108}")
            ax2 = fig.add_axes([0.05, 0.04, 0.90, 0.42])
            ax2.axis('off')
            ax2.text(0.0, 1.0, "\n".join(text_lines),
                     fontsize=8, family='monospace',
                     va='top', ha='left')

            pdf.savefig(fig); plt.close(fig)

    print(f"Wrote: {OUT_PDF}")
    return shape_stats


# ============================================================================
# 6.  WRITE JSON + MARKDOWN
# ============================================================================

def write_text_outputs(shape_stats, n_d9, n_z9, total_multisets):
    n_shapes = len(shape_stats)
    out = {
        "n_z9_cluster_orbits": n_z9,
        "n_d9_groups": n_d9,
        "n_length_spectrum_shapes": n_shapes,
        "n_multisets_total": total_multisets,
        "shapes": [
            {
                "shape_index": i,
                "length_spectrum": list(ss["spectrum"]),
                "spectrum_pretty": spectrum_summary(ss["spectrum"]),
                "n_d9_groups": ss["n_d9_groups"],
                "n_z9_orbits": ss["n_z9_orbits"],
                "n_multisets": ss["n_multisets"],
                "d9_groups": [
                    {
                        "rep_orbit_id": g["rep_orbit_id"],
                        "kind": g["kind"],
                        "z9_orbits_in_group": g["z9_orbits_in_group"],
                        "rep_multiset": g["rep_multiset"],
                    }
                    for g in ss["groups"]
                ],
            }
            for i, ss in enumerate(shape_stats, 1)
        ],
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Wrote: {OUT_JSON}")

    # Markdown
    lines = []
    lines.append("# n=9 cluster: chord-length-spectrum reduction\n")
    lines.append("Coarser than D_9.  Groups multisets by their sorted multiset of "
                 "chord lengths.\n")
    lines.append(f"- Z_9 orbits in cluster: **{n_z9}**")
    lines.append(f"- D_9 groups: **{n_d9}**")
    lines.append(f"- **Length-spectrum shapes: {n_shapes}**")
    lines.append(f"- Total multisets covered: {total_multisets}\n")

    lines.append("## The rule\n")
    lines.append("A chord {a, b} of the 9-gon has LENGTH = min(|a-b|, 9-|a-b|), so "
                 "length ∈ {2, 3, 4}. The LENGTH SPECTRUM of a 6-chord multiset is "
                 "the sorted multiset of chord lengths. Two multisets are "
                 "shape-equivalent iff they have the same length spectrum.\n")
    lines.append("Length is preserved by every rotation and reflection of the "
                 "9-gon, so this is coarser than D_9: every D_9 orbit lives "
                 "inside one length-spectrum shape, but a single shape can "
                 "contain multiple D_9 orbits at different polygon positions.\n")

    lines.append("## All shapes\n")
    lines.append("| Shape # | Length spectrum | # D_9 groups | # Z_9 orbits | # multisets |")
    lines.append("|---:|---|---:|---:|---:|")
    for i, ss in enumerate(shape_stats, 1):
        lines.append(f"| {i} | {spectrum_summary(ss['spectrum'])} | "
                     f"{ss['n_d9_groups']} | {ss['n_z9_orbits']} | "
                     f"{ss['n_multisets']} |")
    lines.append("")

    lines.append("## Shape members (which D_9 groups land in each shape)\n")
    for i, ss in enumerate(shape_stats, 1):
        lines.append(f"### Shape {i} — spectrum {spectrum_summary(ss['spectrum'])}")
        lines.append("")
        lines.append("| D_9 grp rep | Kind | Z_9 orbits | Multiset |")
        lines.append("|---:|---|---|---|")
        for g in ss["groups"]:
            ms = '{' + ', '.join(f"({c[0]},{c[1]})" for c in g["rep_multiset_tuples"]) + '}'
            seed_tag = " ★" if any(o in SEED_ORBITS
                                    for o in g["z9_orbits_in_group"]) else ""
            z9_str = ", ".join(str(o) for o in g["z9_orbits_in_group"])
            lines.append(f"| {g['rep_orbit_id']}{seed_tag} | {g['kind']} | "
                         f"{z9_str} | `{ms}` |")
        lines.append("")
    lines.append("★ = the D_9 group contains a cascade-failure seed orbit "
                 "(22, 46, 88, or 108).")
    with open(OUT_MD, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote: {OUT_MD}")


# ============================================================================
# 7.  MAIN
# ============================================================================

def main():
    d9_data, orbit_size = load_inputs()
    n_d9 = d9_data["n_d9_groups"]
    n_z9 = d9_data["n_z9_cluster_orbits"]

    groups_by_spec = group_by_spectrum(d9_data, orbit_size)
    total_multisets = sum(
        sum(g["total_multisets_in_d9_group"] for g in grps)
        for grps in groups_by_spec.values()
    )

    print(f"D_9 groups: {n_d9}")
    print(f"Distinct length spectra: {len(groups_by_spec)}")
    for spec, grps in sorted(groups_by_spec.items()):
        print(f"  {spectrum_summary(spec):<30}  {len(grps):>2} D_9 groups  "
              f"({sum(len(g['z9_orbits_in_group']) for g in grps)} Z_9 orbits)")

    shape_stats = render_pdf(groups_by_spec, n_d9, n_z9, total_multisets)
    write_text_outputs(shape_stats, n_d9, n_z9, total_multisets)

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
