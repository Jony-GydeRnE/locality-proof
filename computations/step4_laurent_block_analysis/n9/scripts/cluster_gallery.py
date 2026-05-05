"""
Cluster gallery: visualize the 90 n=9 cluster orbit representatives.
================================================================================

Each of the 90 unknowns in the rank-90 cluster linear system corresponds to one
cyclic-orbit representative. This script renders all 90 as filled 9-gons with
red chord overlays — the same visual format as the step-1 survivor gallery —
so each unknown coefficient can be inspected geometrically.

Inputs:
    n9/outputs/orbits_n9.json        — all 113 orbit reps + locality info
    n9/outputs/cluster_partial.json  — list of the 90 cluster orbit IDs

Output:
    n9/outputs/cluster_90_gallery.pdf

Each page shows 6 diagrams (3 rows x 2 cols). The 4 seed orbits that originated
the BFS cluster (22, 46, 88, 108) are flagged in their captions.
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


N = 9
SEED_ORBITS = {22, 46, 88, 108}   # the 4 cascade-failure orbits that seeded the BFS

HERE        = os.path.dirname(os.path.abspath(__file__))
OUTPUTS_DIR = os.path.normpath(os.path.join(HERE, "..", "outputs"))
ORBITS_PATH  = os.path.join(OUTPUTS_DIR, "orbits_n9.json")
CLUSTER_PATH = os.path.join(OUTPUTS_DIR, "cluster_partial.json")
STATUS_PATH  = os.path.join(OUTPUTS_DIR, "n9_locality_status.md")
OUT_PDF      = os.path.join(OUTPUTS_DIR, "cluster_90_gallery.pdf")


def vertex_xy(i, n, R=1.0):
    angle = math.pi / 2 - 2 * math.pi * (i - 1) / n
    return (R * math.cos(angle), R * math.sin(angle))


def parse_locality_tags(path):
    """Pull the Loc tag (CROSSING / DOUBLE_POLE / TRIANGULATION) per orbit_id from the
    consolidated markdown table. Returns dict[int, str]; empty if file missing."""
    tags = {}
    if not os.path.exists(path):
        return tags
    with open(path, "r") as f:
        for line in f:
            if not line.startswith("| ") or "CLUSTER" not in line and "NON_CLUSTER" not in line:
                continue
            parts = [p.strip() for p in line.strip().strip("|").split("|")]
            # columns: orbit, cluster?, Loc, representative, mechanism
            if len(parts) < 3:
                continue
            try:
                oid = int(parts[0])
            except ValueError:
                continue
            tags[oid] = parts[2]
    return tags


def draw_diagram(ax, multiset, n):
    coords = [vertex_xy(i, n) for i in range(1, n + 1)]
    ax.add_patch(MplPolygon(coords, closed=True,
                            facecolor='black', edgecolor='black'))
    chord_counts = {}
    for chord in multiset:
        key = tuple(chord)
        chord_counts[key] = chord_counts.get(key, 0) + 1
    for chord, mult in chord_counts.items():
        i, j = chord
        x1, y1 = vertex_xy(i, n)
        x2, y2 = vertex_xy(j, n)
        for k in range(mult):
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy) or 1.0
            nx, ny = -dy / length, dx / length
            shift = 0.03 * (k - (mult - 1) / 2.0)
            ax.plot([x1 + shift * nx, x2 + shift * nx],
                    [y1 + shift * ny, y2 + shift * ny],
                    color='red', linewidth=2.4)
    # vertex labels (small, white) so geometry is readable
    for v in range(1, n + 1):
        x, y = vertex_xy(v, n, R=1.12)
        ax.text(x, y, str(v), fontsize=6, ha='center', va='center', color='black')
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect('equal')
    ax.axis('off')


def caption_for(orbit_id, multiset, loc_tag):
    ms_str = '{' + ', '.join(f"({c[0]},{c[1]})" for c in multiset) + '}'
    seed_marker = "  [SEED]" if orbit_id in SEED_ORBITS else ""
    loc_str = f"  [{loc_tag}]" if loc_tag else ""
    return f"orbit #{orbit_id}{seed_marker}{loc_str}\n{ms_str}"


def main():
    with open(ORBITS_PATH, "r") as f:
        manifest = json.load(f)
    with open(CLUSTER_PATH, "r") as f:
        cluster = json.load(f)

    cluster_ids = set(cluster["cluster_orbits"])
    if len(cluster_ids) != 90:
        sys.exit(f"Expected 90 cluster orbits, got {len(cluster_ids)}.")

    by_id = {o["orbit_id"]: o for o in manifest["orbits"]}
    loc_tags = parse_locality_tags(STATUS_PATH)

    cluster_orbits = [by_id[oid] for oid in sorted(cluster_ids)]
    print(f"Rendering {len(cluster_orbits)} cluster orbit reps (n={N})...")

    per_page = 6
    cols = 2
    rows = per_page // cols
    n_pages = math.ceil(len(cluster_orbits) / per_page)

    with PdfPages(OUT_PDF) as pdf:
        # cover page
        cover = plt.figure(figsize=(8.5, 11))
        cover.text(0.5, 0.78, "n = 9 cluster gallery",
                   ha='center', fontsize=20, weight='bold')
        cover.text(0.5, 0.72,
                   "All 90 cyclic-orbit representatives whose coefficients are\n"
                   "fixed to zero by the rank-90 cluster linear system\n"
                   "(rows: 506 cousin equations; cols: 90 cluster + 53 external).",
                   ha='center', fontsize=11)
        cover.text(0.5, 0.58,
                   f"4 seed orbits (cascade failures, BFS roots): "
                   + ", ".join(str(s) for s in sorted(SEED_ORBITS)),
                   ha='center', fontsize=10)
        cover.text(0.5, 0.54,
                   "Black 9-gon, red chords = chord multiset of the orbit representative.\n"
                   "Vertex labels 1..9 around the rim. Loc tag = CROSSING or DOUBLE_POLE.",
                   ha='center', fontsize=9, style='italic')
        cover.text(0.5, 0.06,
                   f"Source: n9/outputs/orbits_n9.json + cluster_partial.json\n"
                   f"Pages: {n_pages}, 6 diagrams/page",
                   ha='center', fontsize=8)
        pdf.savefig(cover)
        plt.close(cover)

        for page_idx in range(n_pages):
            start = page_idx * per_page
            page = cluster_orbits[start: start + per_page]
            fig, axes = plt.subplots(rows, cols, figsize=(8.5, 11))
            fig.suptitle(
                f"n=9 cluster orbits  —  page {page_idx + 1} / {n_pages}  "
                f"(orbits {page[0]['orbit_id']}–{page[-1]['orbit_id']})",
                fontsize=11)
            axes_flat = axes.flatten()
            for ax in axes_flat:
                ax.axis('off')
            for ax, orbit in zip(axes_flat, page):
                oid = orbit["orbit_id"]
                ms  = orbit["representative"]
                draw_diagram(ax, ms, N)
                ax.set_title(caption_for(oid, ms, loc_tags.get(oid, "")),
                             fontsize=8, pad=4)
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            pdf.savefig(fig)
            plt.close(fig)
            print(f"  page {page_idx + 1}/{n_pages} done")

    print(f"\nWrote: {OUT_PDF}")
    try:
        if sys.platform == "darwin":
            subprocess.Popen(["open", OUT_PDF])
        elif sys.platform.startswith("linux"):
            subprocess.Popen(["xdg-open", OUT_PDF])
        elif sys.platform == "win32":
            os.startfile(OUT_PDF)  # type: ignore[attr-defined]
        print("Opened PDF in default viewer.")
    except Exception as e:
        print(f"(Could not auto-open PDF: {e})")


if __name__ == "__main__":
    main()
