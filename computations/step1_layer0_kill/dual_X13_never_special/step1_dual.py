"""
================================================================================
DUAL EXPERIMENT  (Item 6 of the master plan)
"What if X_{13} is NEVER the special variable?"
================================================================================

Background you need to understand the question
----------------------------------------------
The Step-1 kill mechanism works zone-by-zone. There are n cyclic zones
Z_{1,3}, Z_{2,4}, ..., Z_{n,2}. At zone Z_{r,r+2} the chord X_{r,r+2} plays
the role of the "special" variable (the one we send to infinity in the
leading-order Laurent argument).

The dual experiment asks:  what if we ALLOW the kill mechanism to use every
zone EXCEPT zone Z_{1,3} ?  Equivalently, what if X_{13} is never allowed to
be the special variable?

Two things can happen:
  1. The other n-1 zones are enough to kill exactly the same multisets.
     Then Z_{1,3} was redundant for Step-1; cyclic symmetry implies any one
     zone is redundant. (This was already observed at n=5 in note 16.)
  2. Some multisets are killed by Z_{1,3} alone and survive when it's
     excluded. These are the "extra survivors" -- the ones that Step 1
     cannot catch without Z_{1,3}, and that Step 2 of the kill mechanism
     (the equality-chain mechanism using the bare relation
        X_{r+1, r-1}  =  -X_{r, r+2}
     i.e., X_{2n} = -X_{13} at Z_{1,3}) must handle instead.

Identifying the extras is therefore the cleanest way to isolate which
non-locals are killed by Step-1 redundancy and which require Step-2's
equality-chain argument. That's the Item 6 deliverable.

What this script does
---------------------
  1. Enumerates every size-(n-3) multiset of chords of the n-gon.
  2. For each multiset:
       (a) checks whether some zone (out of all n) can kill it via Step-1
           -> "baseline survivors";
       (b) checks whether some zone OTHER than the excluded one (default
           r=1, i.e., Z_{1,3}) can kill it
           -> "dual-experiment survivors".
  3. Prints both counts and computes the EXTRAS = (dual survivors) minus
     (baseline survivors). Those extras are precisely the multisets killed
     uniquely by the excluded zone.
  4. Splits each survivor list into triangulations (= the actual tree
     amplitude, expected to survive every zone) versus non-triangulations
     (the "diagrams that refuse to die" -- these are what we care about).
  5. Renders the EXTRA non-triangulation survivors to a PDF, drawn as
     filled n-gons with the chord multiset in red, and auto-opens the PDF.

Usage
-----
    python3 computations/step1_dual.py                # prompts for n
    python3 computations/step1_dual.py 7              # n on the command line
    python3 computations/step1_dual.py 7 --exclude 3  # exclude zone Z_{3,5}
                                                      #   (special chord X_{35})

Output file path:
    computations/visualizing nonlocals that refuse to die/
        nonlocal n=N (dual; Z_R excluded; extras only).pdf
"""

import argparse
import math
import os
import subprocess
import sys
import time
from itertools import combinations_with_replacement

# matplotlib backend "Agg" means "no live window" -- we write a PDF and let
# the OS open it in the default viewer.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Polygon as MplPolygon


# ============================================================================
# 1.  COMBINATORICS OF CHORDS AND ZONES
#     (identical to step1_uncaught.py -- duplicated here for self-contained
#      reading; if you understand the kill mechanism you can skip these.)
# ============================================================================

def normalize(i, j, n):
    """
    Return the canonical name of the chord between vertices i and j of the
    n-gon as a pair (a, b) with a < b in {1, ..., n}.

    Indices may be given mod n (negative, > n, etc.); the function reduces
    them cyclically. This lets us write things like normalize(r+1, r-1, n)
    for the bare chord at zone Z_{r,r+2} without worrying about wrap-around.
    """
    i = ((i - 1) % n) + 1
    j = ((j - 1) % n) + 1
    if i > j:
        i, j = j, i
    return (i, j)


def all_chords(n):
    """
    Every planar chord of the n-gon: pairs (i, j) with 1 <= i < j <= n that
    are NOT polygon edges.

    A pair (i, j) is an edge iff j - i == 1 (consecutive vertices) or
    (i, j) == (1, n) (the wrap-around edge).

    Total chord count = n(n-3)/2.
    """
    return [(i, j) for i in range(1, n + 1)
                   for j in range(i + 1, n + 1)
                   if j - i >= 2 and not (i == 1 and j == n)]


def zone_structure(r, n):
    """
    For zone Z_{r, r+2} return the data the kill criterion needs:

        special  -- the chord X_{r, r+2} that becomes the Laurent variable
                    sent to infinity at this zone.
        bare     -- the chord X_{r+1, r-1} (cyclic) which the 1-zero
                    constraints force to equal -X_{r, r+2}.
        pairs    -- the n-4 (companion, substitute) chord pairs
                       (X_{r, k}, X_{r+1, k}),  k = r+3, ..., r+n-2
                    each linked by   X_{r+1, k}  =  X_{r, k} - X_{r, r+2}.
                    Sending the special to infinity makes EXACTLY ONE side
                    of each pair go to infinity (whichever we choose to);
                    the other stays finite as a free variable.
    """
    special = normalize(r,     r + 2, n)
    bare    = normalize(r + 1, r - 1, n)        # cyclic; non-adjacent chord
    pairs = []
    for offset in range(3, n - 1):              # offsets 3, 4, ..., n-2
        k = ((r - 1 + offset) % n) + 1
        companion  = normalize(r,     k, n)
        substitute = normalize(r + 1, k, n)
        pairs.append((companion, substitute))
    return special, bare, pairs


def killable_at_zone(ms_set, structure):
    """
    Decide whether the multiset (passed in as a Python set for fast
    membership checks) is killable at the zone described by `structure`.

    Two conditions reproduce the leading-order Laurent argument:

      (A) the multiset must NOT contain the special or the bare chord.
          If it did, the corresponding 1/X factor vanishes in the limit
          and the monomial drops out at leading order, producing NO
          equation about that coefficient -- so it cannot be killed at
          this zone via Step-1.

      (B) for every (companion, substitute) pair, the multiset must NOT
          contain BOTH chords of the pair.
          If it did, no consistent sub-limit holds either side of the
          pair finite while sending the special to infinity, so again the
          monomial cannot be isolated at leading order.

    If both (A) and (B) hold, we can pick a sub-limit consistent with the
    multiset and the linear-independence-of-monomials argument kills the
    coefficient.
    """
    special, bare, pairs = structure
    if special in ms_set or bare in ms_set:     # condition (A)
        return False
    for companion, substitute in pairs:         # condition (B)
        if companion in ms_set and substitute in ms_set:
            return False
    return True


# ============================================================================
# 2.  GEOMETRY -- CROSSING TEST AND TRIANGULATION TEST
# ============================================================================

def chords_cross(c1, c2):
    """
    Two chords (a,b) and (c,d) of the convex n-gon (vertices in cyclic
    order 1..n) cross iff their four endpoints alternate around the
    boundary -- equivalently exactly one endpoint of (c,d) lies strictly
    between a and b in the linear order a < ... < b.
    Sharing an endpoint => no crossing (chords meet only at a vertex).
    """
    a, b = c1
    c, d = c2
    if {a, b} & {c, d}:
        return False
    inside = lambda v: a < v < b
    return inside(c) != inside(d)


def is_triangulation(ms):
    """
    A size-(n-3) multiset is a triangulation of the n-gon iff
        (i)  all its chords are distinct (no double poles), and
        (ii) no two chords cross.
    Then the chords carve the polygon into n-2 triangles -- this is what
    we call a tree-amplitude diagram. There are Catalan-many of them.
    """
    if len(set(ms)) != len(ms):
        return False
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            if chords_cross(ms[i], ms[j]):
                return False
    return True


# ============================================================================
# 3.  DIAGRAMS  (same renderer as step1_uncaught.py)
# ============================================================================

def vertex_xy(i, n, R=1.0):
    """Cartesian position of vertex i of a regular n-gon inscribed in a
    circle of radius R (vertex 1 at top, going clockwise)."""
    angle = math.pi / 2 - 2 * math.pi * (i - 1) / n
    return (R * math.cos(angle), R * math.sin(angle))


def draw_diagram(ax, multiset, n):
    """Filled black n-gon with red chords for the multiset.  Multi-occurrence
    chords (double poles) are shown as parallel offset copies."""
    coords = [vertex_xy(i, n) for i in range(1, n + 1)]
    ax.add_patch(MplPolygon(coords, closed=True,
                            facecolor='black', edgecolor='black'))
    chord_counts = {}
    for chord in multiset:
        chord_counts[chord] = chord_counts.get(chord, 0) + 1
    for chord, mult in chord_counts.items():
        i, j = chord
        x1, y1 = vertex_xy(i, n)
        x2, y2 = vertex_xy(j, n)
        for k in range(mult):
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy) or 1.0
            nx, ny = -dy / length, dx / length      # unit normal
            shift = 0.025 * (k - (mult - 1) / 2.0)   # symmetric offsets
            ax.plot([x1 + shift * nx, x2 + shift * nx],
                    [y1 + shift * ny, y2 + shift * ny],
                    color='red', linewidth=2.4)
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_aspect('equal')
    ax.axis('off')
    label_str = '{' + ', '.join('{' + f'{c[0]},{c[1]}' + '}' for c in multiset) + '}'
    ax.text(1.4, 0, f"$a_{{{label_str}}}$", fontsize=7, va='center', ha='left')


def build_figures(multisets, n, title=None, cols=2, per_page=8):
    """Pack diagrams onto letter-size pages, two columns per row.
    Returns a list of figures the caller dumps into a PDF."""
    if not multisets:
        return []
    rows = math.ceil(per_page / cols)
    figs = []
    for page_idx, page_start in enumerate(range(0, len(multisets), per_page)):
        page = multisets[page_start: page_start + per_page]
        fig, axes = plt.subplots(rows, cols, figsize=(8.5, 11))
        if title:
            suffix = "" if page_idx == 0 else f"  (page {page_idx + 1})"
            fig.suptitle(title + suffix, fontsize=11)
        axes_flat = axes.flatten() if rows * cols > 1 else [axes]
        for ax in axes_flat:
            ax.axis('off')
        for ax, ms in zip(axes_flat, page):
            draw_diagram(ax, ms, n)
        fig.tight_layout(rect=[0, 0, 0.95, 0.96])
        figs.append(fig)
    return figs


# ============================================================================
# 4.  THE DUAL EXPERIMENT
#     -- enumerate, classify, compare, render
# ============================================================================

def zone_label(r, n):
    """Pretty 'Z_{r,s}' label string for zone Z_{r, r+2} (cyclic-aware)."""
    a, b = normalize(r, r + 2, n)
    return f"Z_{{{a},{b}}}"


def parse_args():
    """Get n and the excluded zone index from CLI; prompt for n if missing."""
    parser = argparse.ArgumentParser(
        description="Step-1 kill mechanism with one zone excluded."
    )
    parser.add_argument("n", nargs="?", type=int, default=None,
                        help="Number of legs of the polygon (>= 4).")
    parser.add_argument("--exclude", type=int, default=1,
                        help="Zone index r to exclude. Default is 1, "
                             "which is Z_{1,3} (the zone where X_{13} is "
                             "the special variable).")
    args = parser.parse_args()
    if args.n is None:
        args.n = int(input("Enter n (>= 4, recommended <= 9): ").strip())
    return args.n, args.exclude


def run_dual(n, excluded_zone):
    """
    Main routine:
      - enumerate every size-(n-3) multiset of chords
      - classify each as killable (or not) under (a) the full-zone Step-1
        (baseline) and (b) the all-zones-except-excluded Step-1 (dual).
      - compute and report the EXTRAS = dual survivors that the baseline
        does NOT have, i.e., the multisets uniquely killed by the
        excluded zone.
      - render the non-triangulation extras to a PDF and auto-open it.
    """
    N = n - 3                                   # multiset size
    chords     = all_chords(n)
    structures = [zone_structure(r, n) for r in range(1, n + 1)]

    # the structure of the EXCLUDED zone (so we can print which special
    # chord it corresponds to, for human readability)
    sp_excl, _, _ = structures[excluded_zone - 1]

    print(f"\nn = {n}, multiset size N = {N}, |chords| = {len(chords)}")
    print(f"Excluding zone {zone_label(excluded_zone, n)} "
          f"(its special chord is X_{{{sp_excl[0]}{sp_excl[1]}}}).")
    print(f"Step-1 kill is run with the remaining {n-1} zones.\n")

    # ----- enumerate ----------------------------------------------------
    print("Enumerating multisets...")
    t0 = time.time()
    survivors_full = []   # uncaught with ALL n zones (baseline)
    survivors_dual = []   # uncaught with n-1 zones (dual experiment)
    for ms in combinations_with_replacement(chords, N):
        ms_set = set(ms)
        # baseline: any zone can kill it?
        killed_full = any(killable_at_zone(ms_set, s) for s in structures)
        # dual: any zone EXCEPT the excluded one can kill it?
        killed_dual = any(killable_at_zone(ms_set, structures[r - 1])
                          for r in range(1, n + 1) if r != excluded_zone)
        if not killed_full:
            survivors_full.append(ms)
        if not killed_dual:
            survivors_dual.append(ms)
    dt = time.time() - t0
    print(f"  ({dt:.2f} s)\n")

    # ----- classify -----------------------------------------------------
    tri_full    = [m for m in survivors_full if     is_triangulation(m)]
    nontri_full = [m for m in survivors_full if not is_triangulation(m)]
    tri_dual    = [m for m in survivors_dual if     is_triangulation(m)]
    nontri_dual = [m for m in survivors_dual if not is_triangulation(m)]

    # extras = multisets the baseline catches but the dual misses
    full_set = set(survivors_full)              # set for O(1) membership
    extras       = [m for m in survivors_dual if m not in full_set]
    extras_nontri = [m for m in extras if not is_triangulation(m)]

    # ----- report -------------------------------------------------------
    print(f"Step-1 baseline (all {n} zones):")
    print(f"   {len(survivors_full):6d} survivors  "
          f"=  {len(tri_full):4d} triangulations  +  "
          f"{len(nontri_full):4d} non-triangulations")
    print(f"Step-1 dual  ({n-1} zones; "
          f"{zone_label(excluded_zone, n)} excluded):")
    print(f"   {len(survivors_dual):6d} survivors  "
          f"=  {len(tri_dual):4d} triangulations  +  "
          f"{len(nontri_dual):4d} non-triangulations")
    print()
    print(f"EXTRAS  (uniquely killed by {zone_label(excluded_zone, n)}):  "
          f"{len(extras)}  total")
    print(f"  -- of which non-triangulations:  {len(extras_nontri)}")
    print()
    print("INTERPRETATION")
    print("  The extras are the non-locals that the OTHER zones cannot")
    print("  catch via Step 1.  In the master plan they are the multisets")
    print("  that Step 2 of the kill mechanism (the bare-relation cyclic-")
    print("  shift / equality-chain argument) must handle when the")
    print("  excluded zone is taken off the table.")

    if extras_nontri and len(extras_nontri) <= 30:
        print(f"\nExtra non-triangulation survivors ({len(extras_nontri)}):")
        for m in extras_nontri:
            print(f"  {m}")
    elif extras_nontri:
        print(f"\n(omitting print of {len(extras_nontri)} multisets; see PDF)")

    # ----- render PDF ---------------------------------------------------
    here    = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(here, "outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir,
        f"nonlocal n={n} (dual; {zone_label(excluded_zone, n)} excluded; extras only).pdf"
    )

    print(f"\nRendering {len(extras_nontri)} diagram(s)...")
    figs = build_figures(extras_nontri, n,
        title=f"n={n}: non-tri multisets uniquely killed by "
              f"{zone_label(excluded_zone, n)} "
              f"(special X_{{{sp_excl[0]}{sp_excl[1]}}})")
    with PdfPages(out_path) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Wrote: {out_path}")

    # auto-open
    try:
        if sys.platform == "darwin":
            subprocess.Popen(["open", out_path])
        elif sys.platform.startswith("linux"):
            subprocess.Popen(["xdg-open", out_path])
        elif sys.platform == "win32":
            os.startfile(out_path)              # type: ignore[attr-defined]
        print("Opened PDF in default viewer.")
    except Exception as e:
        print(f"(Could not auto-open PDF: {e})")


def main():
    n, excluded = parse_args()
    if n < 4:
        sys.exit(f"n = {n} too small; need n >= 4.")
    if not 1 <= excluded <= n:
        sys.exit(f"--exclude must be in [1, {n}].")
    run_dual(n, excluded)


if __name__ == "__main__":
    main()
