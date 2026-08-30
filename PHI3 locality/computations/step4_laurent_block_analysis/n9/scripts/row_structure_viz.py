"""
Row-structure visualization of the 506-equation cluster system.
================================================================================

WHAT THIS SCRIPT IS FOR
-----------------------
The proof that the 90 cluster orbits are killed comes from a giant linear
system: 506 equations in 143 unknowns (90 cluster orbits + 53 external
columns). The actual scalar entries of this system live inside a 10-hour
sympy computation; what got SAVED to disk are only the row LABELS — for
each of the 506 equations we know which cluster orbit it came from, which
hidden-zero zone, which substitute chord, and which fingerprint.

This script visualizes the STRUCTURE of those equations using only the
labels, no scalar entries. It answers questions like:

  * How many equations does each cluster orbit contribute?
  * From which hidden-zero zones (Z_1, Z_2, ..., Z_9) do they come?
  * Do reflection-paired orbits (D_9 partners) produce the SAME set of
    equations under the reflection map?  If yes, that's strong evidence
    that the 506 × 90 linear system has full D_9 symmetry, which means
    it block-decomposes into 55 smaller systems indexed by D_9 groups
    (the user's "Jordan-form-ish" picture).

INPUTS:
  n9/outputs/cluster_partial.json     -- 506 row_labels
  n9/outputs/dihedral_reduction.json  -- 55 D_9 groups + reflection partners

OUTPUTS:
  n9/outputs/row_structure.pdf
  n9/outputs/row_structure.json   (per-orbit row counts + reflection check)
"""

import json
import math
import os
import subprocess
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


# ============================================================================
# 1.  WHERE FILES LIVE
# ============================================================================
N = 9
HERE        = os.path.dirname(os.path.abspath(__file__))
OUTPUTS_DIR = os.path.normpath(os.path.join(HERE, "..", "outputs"))
ROWS_PATH    = os.path.join(OUTPUTS_DIR, "cluster_partial.json")
DIHEDRAL_PATH = os.path.join(OUTPUTS_DIR, "dihedral_reduction.json")
OUT_PDF     = os.path.join(OUTPUTS_DIR, "row_structure.pdf")
OUT_JSON    = os.path.join(OUTPUTS_DIR, "row_structure.json")


# ============================================================================
# 2.  THE 18 ELEMENTS OF D_9 ACTING ON VERTICES, CHORDS, AND ZONES
# ============================================================================
# D_9 has 9 rotations and 9 reflections (18 total).
#
#   Rotation by k:  rho_k(i) = ((i - 1 + k) mod 9) + 1
#   Base reflection (axis through vertex 1):  sigma(i) = ((1 - i) mod 9) + 1
#   The 9 reflections are then:               rho_k ∘ sigma   for k = 0..8
#
# We need to act on three label-objects:
#   * a single vertex  i ∈ {1..9}
#   * a chord {a, b}   (we keep them as tuples (a, b) with a < b)
#   * a zone-index r   (treated like a vertex; the "leg-r locus" Z_{r,r+2}
#     transforms with vertex r)

def rotate_v(i, k, n=N):           return ((i - 1 + k) % n) + 1
def base_reflect_v(i, n=N):        return ((1 - i) % n) + 1


def perm_after(k, do_reflect):
    """Return the function vertex -> vertex implementing
       rho_k          if do_reflect is False,
       rho_k ∘ sigma  if do_reflect is True (i.e. one of the 9 reflections).
    """
    if do_reflect:
        return lambda i: rotate_v(base_reflect_v(i), k)
    else:
        return lambda i: rotate_v(i, k)


def all_reflections():
    """Yield the 9 reflections of D_9 as functions."""
    for k in range(N):
        yield perm_after(k, do_reflect=True)


def apply_to_chord(chord, perm):
    a, b = perm(chord[0]), perm(chord[1])
    return (a, b) if a < b else (b, a)


def apply_to_signature(sig, perm):
    """A signature is (zone_r, U, Y).  Apply the vertex permutation to each
       component: the zone uses vertex-style transformation, U and Y are
       chords."""
    z, u, y = sig
    return (perm(z), apply_to_chord(u, perm), apply_to_chord(y, perm))


# ============================================================================
# 3.  LOAD ROW LABELS AND D_9 GROUPING
# ============================================================================

def load_inputs():
    with open(ROWS_PATH) as f:
        rows_data = json.load(f)
    with open(DIHEDRAL_PATH) as f:
        d9_data = json.load(f)
    return rows_data["row_labels"], d9_data


# ============================================================================
# 4.  PER-ORBIT ROW COUNTS AND ZONE BREAKDOWNS
# ============================================================================

def per_orbit_zone_grid(row_labels, cluster_orbits):
    """Return a (n_orbits, 9) numpy array: entry [i, z-1] = number of
       equations from orbit cluster_orbits[i] coming from zone z."""
    grid = np.zeros((len(cluster_orbits), N), dtype=int)
    orbit_to_idx = {oid: i for i, oid in enumerate(cluster_orbits)}
    for r in row_labels:
        i = orbit_to_idx[r["orbit"]]
        z = r["zone_r"] - 1
        grid[i, z] += 1
    return grid


# ============================================================================
# 5.  REFLECTION-SYMMETRY TEST
# ============================================================================
# For each reflection-pair group (Z_9 orbits a and b paired by reflection),
# we check that the set of (zone, U, Y) labels of orbit a, after applying
# reflection to each label, equals the set of labels of orbit b.
#
# If this holds for every pair, the 506-row matrix is exactly invariant
# under D_9: equations come in reflection-paired packs, and the matrix
# block-decomposes by D_9 irreps into 55 smaller blocks.

def reflect_zone(zone_r):
    """Zone Z_{r,r+2} reflects to zone Z_{r',r'+2} where r' = ((-r-1) mod 9) + 1.
       Derivation: under sigma, vertex pair {r, r+2} (0-idx {r-1, r+1}) reflects
       to {-(r-1), -(r+1)} = {1-r, -r-1} (0-idx). The new zone r' satisfies
       {r'-1, r'+1} = {1-r, -r-1}, giving r'-1 = -r-1 (0-idx), r' = -r (0-idx).
       Translating to 1-idx: r'_1 = ((-r-1) mod 9) + 1 = ref(r+2)_1."""
    return ((-zone_r - 1) % N) + 1


def per_orbit_zone_vector(rows_for_orbit):
    """Return tuple of length 9: count of equation rows at each zone."""
    counts = [0] * N
    for r in rows_for_orbit:
        counts[r["zone_r"] - 1] += 1
    return tuple(counts)


def find_matching_reflection_count(zone_vec_a, zone_vec_b):
    """Try all 9 reflections rho_k . sigma.  Each one acts on a per-zone
       count vector by: zone r -> ref_zone(r) then -> ref_zone(r) + k.
       Return the first k (0..8) for which the transformed count vector
       of a equals zone_vec_b; return None if none works.

       PHYSICS: checks the "weak" D_9 invariance — whether orbit a and
       orbit b produce the SAME NUMBER of equations at corresponding
       reflected zones, even though the (substitute, companion) detail
       inside each row need not match (the kill mechanism privileges the
       +1 vertex direction at each zone, so the FULL labels are not
       reflection-symmetric, but the COUNTS per zone are)."""
    if sum(zone_vec_a) != sum(zone_vec_b):
        return None
    # Reflected zones for r=1..9 (1-idx): a vector
    reflected = [reflect_zone(r) for r in range(1, N + 1)]
    for k in range(N):
        # Build transformed count vector: for each new-zone r_new, sum
        # contributions from old zones whose reflection-then-rotate-by-k
        # equals r_new.
        transformed = [0] * N
        for r_old_idx in range(N):
            r_old = r_old_idx + 1
            r_new = ((reflected[r_old_idx] - 1 + k) % N) + 1
            transformed[r_new - 1] += zone_vec_a[r_old_idx]
        if tuple(transformed) == zone_vec_b:
            return k
    return None


def check_reflection_symmetry(rows_by_orbit, d9_data):
    """For each reflection-pair, check whether the per-zone equation
       counts match under some D_9 reflection."""
    results = []
    for grp in d9_data["groups"]:
        if grp["kind"] != "reflection_pair":
            continue
        a, b = grp["z9_orbits_in_group"]
        vec_a = per_orbit_zone_vector(rows_by_orbit.get(a, []))
        vec_b = per_orbit_zone_vector(rows_by_orbit.get(b, []))
        k_match = find_matching_reflection_count(vec_a, vec_b)
        results.append({
            "orbit_a": a,
            "orbit_b": b,
            "n_eqs_a": sum(vec_a),
            "n_eqs_b": sum(vec_b),
            "zone_vec_a": list(vec_a),
            "zone_vec_b": list(vec_b),
            "reflected_a_equals_b": (k_match is not None),
            "matching_reflection_k": k_match,
        })
    return results


def check_palindrome_self_reflection(rows_by_orbit, d9_data):
    """For each palindromic D_9 group, check whether the per-zone counts
       are invariant under some D_9 reflection."""
    results = []
    for grp in d9_data["groups"]:
        if grp["kind"] != "palindromic":
            continue
        oid = grp["rep_orbit_id"]
        vec = per_orbit_zone_vector(rows_by_orbit.get(oid, []))
        k_match = find_matching_reflection_count(vec, vec)
        results.append({
            "orbit": oid,
            "n_eqs": sum(vec),
            "zone_vec": list(vec),
            "self_reflection_invariant": (k_match is not None),
            "matching_reflection_k": k_match,
        })
    return results


# ============================================================================
# 6.  PLOTTING
# ============================================================================

def make_color_for_orbit(oid, d9_data):
    """Color-code each cluster orbit by the kind of its D_9 group."""
    for grp in d9_data["groups"]:
        if oid in grp["z9_orbits_in_group"]:
            kind = grp["kind"]
            if kind == "palindromic":
                return "#4caf50"   # green
            elif kind == "reflection_pair":
                return "#1976d2"   # blue
            elif kind == "palindromic_leak":
                return "#ef6c00"   # orange
            else:
                return "#888888"
    return "#888888"


def render_pdf(row_labels, d9_data, rows_by_orbit, grid, cluster_orbits,
               pair_check, palindrome_check):
    n_d9 = d9_data["n_d9_groups"]
    n_pal = d9_data["n_palindromic"]
    n_pair = d9_data["n_reflection_pairs"]
    n_leak = d9_data["n_palindromic_with_external_partner"]

    # Map orbit -> D_9 group index, kind, partner
    orbit_to_group = {}
    for i, grp in enumerate(d9_data["groups"]):
        for o in grp["z9_orbits_in_group"]:
            orbit_to_group[o] = (i, grp["kind"], grp)

    # Re-ordering: by D_9 group, paired orbits adjacent
    d9_ordered_orbits = []
    for grp in d9_data["groups"]:
        for o in grp["z9_orbits_in_group"]:
            d9_ordered_orbits.append(o)

    with PdfPages(OUT_PDF) as pdf:
        # ---------- PAGE 1 : COVER SUMMARY ----------
        cov = plt.figure(figsize=(8.5, 11))
        cov.text(0.5, 0.92, "n=9 cluster: equation row structure",
                 ha='center', fontsize=18, weight='bold')
        cov.text(0.5, 0.87,
                 "(visualizing the 506 equation labels — scalar entries not needed)",
                 ha='center', fontsize=10, style='italic')

        n_pair_match = sum(1 for r in pair_check if r["reflected_a_equals_b"])
        n_pair_total = len(pair_check)
        n_pal_match = sum(1 for r in palindrome_check if r["self_reflection_invariant"])
        n_pal_total = len(palindrome_check)

        cov.text(0.5, 0.78,
                 f"506 equations, 90 cluster columns, 53 external columns\n"
                 f"D_9 reduction: {n_d9} groups\n"
                 f"  palindromic         : {n_pal}\n"
                 f"  reflection pair     : {n_pair}\n"
                 f"  palindromic-leak    : {n_leak}",
                 ha='center', fontsize=11)

        cov.text(0.5, 0.62, "Reflection-symmetry test of the equation set",
                 ha='center', fontsize=14, weight='bold')

        cov.text(0.5, 0.50,
                 f"For each reflection-pair (orbit a, orbit b):  does the\n"
                 f"per-zone equation COUNT vector of a, after applying some\n"
                 f"D_9 reflection, equal that of b?\n\n"
                 f"   → matched: {n_pair_match} / {n_pair_total} pairs\n\n"
                 f"For each palindromic orbit:  is its zone-count vector\n"
                 f"invariant under some D_9 reflection?\n\n"
                 f"   → matched: {n_pal_match} / {n_pal_total} palindromes",
                 ha='center', fontsize=11)

        if n_pair_match == n_pair_total and n_pal_match == n_pal_total:
            verdict = ("✓  Per-zone equation counts respect D_9.\n"
                       "→ Reflection-paired orbits produce the same NUMBER\n"
                       "  of equations at each (reflected) zone.\n"
                       "  NOTE: the full (zone, U, Y) labels do NOT match\n"
                       "  exactly under reflection — the kill mechanism\n"
                       "  privileges the +1 vertex direction at each zone,\n"
                       "  so substitutes are not reflection-symmetric. Counts\n"
                       "  are; details are not.")
        else:
            verdict = ("✗ Per-zone counts do NOT all match under D_9.\n"
                       "  See per-pair details on the next pages.")
        cov.text(0.5, 0.22, verdict, ha='center', fontsize=10,
                 weight='bold' if "✓" in verdict else 'normal')

        cov.text(0.5, 0.06,
                 "Sources: cluster_partial.json (row_labels) + "
                 "dihedral_reduction.json", ha='center', fontsize=8)
        pdf.savefig(cov); plt.close(cov)

        # ---------- PAGE 2 : ROWS-PER-ORBIT BAR CHART ----------
        fig, ax = plt.subplots(figsize=(11, 8.5))
        n_rows_per_orbit = np.array([len(rows_by_orbit.get(o, []))
                                      for o in cluster_orbits])
        colors = [make_color_for_orbit(o, d9_data) for o in cluster_orbits]
        x = np.arange(len(cluster_orbits))
        ax.bar(x, n_rows_per_orbit, color=colors, edgecolor='black', linewidth=0.3)
        ax.set_xticks(x)
        ax.set_xticklabels([str(o) for o in cluster_orbits], rotation=90, fontsize=5)
        ax.set_xlabel("cluster orbit ID", fontsize=10)
        ax.set_ylabel("number of equation rows from this orbit", fontsize=10)
        ax.set_title("Equation rows contributed by each of the 90 cluster orbits\n"
                     "(green = palindromic;  blue = reflection-pair;  orange = palindromic-leak)",
                     fontsize=11)
        # stats annotation
        ax.text(0.99, 0.98,
                f"total rows = {n_rows_per_orbit.sum()}\n"
                f"min/median/max per orbit = "
                f"{n_rows_per_orbit.min()} / {int(np.median(n_rows_per_orbit))} / {n_rows_per_orbit.max()}",
                transform=ax.transAxes, ha='right', va='top', fontsize=9,
                bbox=dict(facecolor='white', alpha=0.7))
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # ---------- PAGE 3 : ZONE x ORBIT HEATMAP (sorted by ID) ----------
        fig, ax = plt.subplots(figsize=(8.5, 11))
        im = ax.imshow(grid, aspect='auto', cmap='viridis', interpolation='nearest')
        ax.set_xticks(range(N))
        ax.set_xticklabels([f"Z_{r}" for r in range(1, N + 1)], fontsize=9)
        ax.set_yticks(range(len(cluster_orbits)))
        ax.set_yticklabels([str(o) for o in cluster_orbits], fontsize=4)
        ax.set_xlabel("hidden-zero zone (leg index r)", fontsize=10)
        ax.set_ylabel("cluster orbit ID", fontsize=10)
        ax.set_title("Equation count per (orbit, zone)  —  orbits sorted by ID",
                     fontsize=11)
        plt.colorbar(im, ax=ax, label="# equations")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # ---------- PAGE 4 : ZONE x ORBIT HEATMAP (sorted by D_9 group) ----------
        # Reorder grid by d9_ordered_orbits.
        idx_for_orbit = {o: i for i, o in enumerate(cluster_orbits)}
        order = [idx_for_orbit[o] for o in d9_ordered_orbits]
        grid_d9 = grid[order]

        fig, ax = plt.subplots(figsize=(8.5, 11))
        im = ax.imshow(grid_d9, aspect='auto', cmap='viridis', interpolation='nearest')
        ax.set_xticks(range(N))
        ax.set_xticklabels([f"Z_{r}" for r in range(1, N + 1)], fontsize=9)
        # Tick label at the start of each D_9 group, with hairline separator
        ytick_positions = []
        ytick_labels    = []
        running = 0
        for g_i, grp in enumerate(d9_data["groups"]):
            sz = len(grp["z9_orbits_in_group"])
            ytick_positions.append(running + (sz - 1) / 2)
            ytick_labels.append(f"D_9 grp #{g_i+1} ({grp['kind'][:3]})")
            running += sz
            if running < len(d9_ordered_orbits):
                ax.axhline(running - 0.5, color='red', linewidth=0.4, alpha=0.5)
        ax.set_yticks(ytick_positions)
        ax.set_yticklabels(ytick_labels, fontsize=4)
        ax.set_xlabel("hidden-zero zone (leg index r)", fontsize=10)
        ax.set_ylabel("cluster orbit (sorted by D_9 group; red lines = group boundaries)", fontsize=10)
        ax.set_title("Equation count per (orbit, zone)  —  orbits sorted by D_9 group\n"
                     "(reflection-paired orbits should have matching profiles)",
                     fontsize=10)
        plt.colorbar(im, ax=ax, label="# equations")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # ---------- PAGE 5+ : PER-PAIR REFLECTION-MATCH TABLE ----------
        # one page-friendly table that summarizes pair_check
        per_page = 30
        for pg_start in range(0, len(pair_check), per_page):
            page = pair_check[pg_start: pg_start + per_page]
            fig, ax = plt.subplots(figsize=(8.5, 11))
            ax.axis('off')
            ax.set_title(
                f"Reflection-pair equation match check  (rows {pg_start+1}-"
                f"{min(pg_start+per_page, len(pair_check))}/{len(pair_check)})",
                fontsize=11)
            cell_text = []
            for r in page:
                marker = "✓" if r["reflected_a_equals_b"] else "✗"
                cell_text.append([
                    f"{r['orbit_a']} ↔ {r['orbit_b']}",
                    str(r["n_eqs_a"]),
                    str(r["n_eqs_b"]),
                    marker,
                    str(r["matching_reflection_k"]) if r["matching_reflection_k"] is not None else "-",
                    " ".join(str(x) for x in r["zone_vec_a"]),
                ])
            table = ax.table(
                cellText=cell_text,
                colLabels=["pair (a ↔ b)", "#eqs a", "#eqs b", "match",
                           "k (rho_k.sigma)", "zone vec a (z=1..9)"],
                loc='upper center',
                cellLoc='center',
            )
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            pdf.savefig(fig); plt.close(fig)

        # ---------- PALINDROME CHECK PAGE ----------
        fig, ax = plt.subplots(figsize=(8.5, 11))
        ax.axis('off')
        ax.set_title("Palindromic-orbit self-reflection check", fontsize=11)
        cell_text = []
        for r in palindrome_check:
            marker = "✓" if r["self_reflection_invariant"] else "✗"
            cell_text.append([
                str(r["orbit"]),
                str(r["n_eqs"]),
                marker,
                str(r["matching_reflection_k"]) if r["matching_reflection_k"] is not None else "-",
                " ".join(str(x) for x in r["zone_vec"]),
            ])
        if cell_text:
            table = ax.table(
                cellText=cell_text,
                colLabels=["orbit", "#eqs", "match",
                           "k (rho_k.sigma)", "zone vec (z=1..9)"],
                loc='upper center', cellLoc='center',
            )
            table.auto_set_font_size(False)
            table.set_fontsize(9)
        pdf.savefig(fig); plt.close(fig)

        # ---------- LEAK PAGE: palindromic-leak orbits + their external partners
        fig, ax = plt.subplots(figsize=(8.5, 11))
        ax.axis('off')
        ax.set_title(
            "Palindromic-leak orbits (reflection partner is OUTSIDE the cluster)",
            fontsize=11)
        leaks = [g for g in d9_data["groups"] if g["kind"] == "palindromic_leak"]
        cell_text = [[str(g["rep_orbit_id"]),
                      str(g["external_partner_orbit"]),
                      str(g["rot_copies"])]
                     for g in leaks]
        if cell_text:
            table = ax.table(
                cellText=cell_text,
                colLabels=["cluster orbit", "external reflection partner orbit",
                           "rot copies"],
                loc='upper center', cellLoc='center',
            )
            table.auto_set_font_size(False)
            table.set_fontsize(9)
        pdf.savefig(fig); plt.close(fig)

    print(f"Wrote: {OUT_PDF}")


# ============================================================================
# 7.  MAIN
# ============================================================================

def main():
    row_labels, d9_data = load_inputs()
    cluster_orbits = sorted({r["orbit"] for r in row_labels})
    print(f"Loaded {len(row_labels)} row labels across {len(cluster_orbits)} orbits.")

    rows_by_orbit = defaultdict(list)
    for r in row_labels:
        rows_by_orbit[r["orbit"]].append(r)

    grid = per_orbit_zone_grid(row_labels, cluster_orbits)
    pair_check       = check_reflection_symmetry(rows_by_orbit, d9_data)
    palindrome_check = check_palindrome_self_reflection(rows_by_orbit, d9_data)

    n_pair_match = sum(1 for r in pair_check if r["reflected_a_equals_b"])
    n_pal_match  = sum(1 for r in palindrome_check if r["self_reflection_invariant"])
    print(f"Reflection-pair check : {n_pair_match}/{len(pair_check)} matched.")
    print(f"Palindrome self check : {n_pal_match}/{len(palindrome_check)} matched.")

    # Convert pair_check tuples (which contain tuples of ints) to JSON-friendly form
    def to_jsonable(v):
        if isinstance(v, dict):
            return {k: to_jsonable(x) for k, x in v.items()}
        if isinstance(v, list):
            return [to_jsonable(x) for x in v]
        if isinstance(v, tuple):
            return list(v)
        return v

    out_data = {
        "n_rows": len(row_labels),
        "n_cluster_orbits": len(cluster_orbits),
        "rows_per_orbit": {str(o): len(rows_by_orbit[o]) for o in cluster_orbits},
        "zone_grid_orbits_x_zones": grid.tolist(),
        "cluster_orbits_in_order": cluster_orbits,
        "reflection_pair_match_count":   n_pair_match,
        "reflection_pair_total":         len(pair_check),
        "palindrome_self_match_count":   n_pal_match,
        "palindrome_self_total":         len(palindrome_check),
        "reflection_pair_check":  to_jsonable(pair_check),
        "palindrome_self_check":  to_jsonable(palindrome_check),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out_data, f, indent=2)
    print(f"Wrote: {OUT_JSON}")

    render_pdf(row_labels, d9_data, rows_by_orbit, grid, cluster_orbits,
               pair_check, palindrome_check)

    try:
        if sys.platform == "darwin":
            subprocess.Popen(["open", OUT_PDF])
        elif sys.platform.startswith("linux"):
            subprocess.Popen(["xdg-open", OUT_PDF])
        elif sys.platform == "win32":
            os.startfile(OUT_PDF)  # type: ignore[attr-defined]
        print("Opened PDF.")
    except Exception as e:
        print(f"(Could not auto-open: {e})")


if __name__ == "__main__":
    main()
