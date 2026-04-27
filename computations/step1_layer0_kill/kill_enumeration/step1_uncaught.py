"""
Step-1 (Layer-0) kill mechanism enumeration at any n, with diagram output.
================================================================================

Usage:
    python3 computations/step1_uncaught.py              # prompts for n
    python3 computations/step1_uncaught.py 8            # n given on the CLI

WHAT THE PROGRAM ANSWERS
------------------------
We consider the most general rational ansatz B with mass-dimension -2(n-3)
and only planar simple poles, expanded as a sum over multisets of n-3 chords:

      B  =  sum over (chord_1, ..., chord_{n-3})   a_{chord_1 ... chord_{n-3}}
                                                  -----------------------------
                                                    X_{chord_1} ... X_{chord_{n-3}}

The "Step-1" (a.k.a. Layer-0, leading-order) kill mechanism takes a single
hidden-zero locus Z_{r, r+2} (the leg-r locus, special variable X_{r,r+2}),
sends X_{r,r+2} -> infinity, and reads off which coefficients a_{...} must
vanish from the leading-order survival argument:

   - the special X_{r,r+2} and its "bare" companion X_{r+1, r-1} go to
     +/- infinity, so any monomial whose multiset contains either of them
     drops out of the leading order entirely (no equation produced);
   - for every (companion, substitute) pair (X_{r,k}, X_{r+1,k}) with
     X_{r+1,k} = X_{r,k} - X_{r,r+2}, we may freely choose to hold either
     side arbitrary; the OTHER side then goes to infinity. So a monomial
     using BOTH members of a pair has no consistent free-variable limit
     and is also dropped;
   - the remaining "good" monomials -- those entirely composed of free
     chords (the inner-triangle F^△ and one chord from each pair) -- form
     a sum of distinct rational monomials in algebraically independent
     free variables. Linear independence kills each coefficient.

So at each zone, a multiset is KILLABLE (Step-1) iff
   (A) it contains neither the special nor the bare;
   (B) for every (companion, substitute) pair, it does not contain BOTH.

A multiset is UNCAUGHT by Step-1 iff no zone Z_{r, r+2} (r = 1, ..., n)
makes it killable. Among the uncaught multisets, the triangulations of the
n-gon are EXPECTED to survive (they are the actual amplitude). Any
NON-triangulation uncaught multiset is a coefficient that the Step-1
mechanism alone cannot eliminate -- those are the diagrams "that refuse to
die" and must be killed by a higher-layer mechanism (Step-2, equality
chains, etc.).

OUTPUT
------
- Prints a summary: total multisets, uncaught count, triangulation vs
  non-triangulation breakdown.
- Writes a PDF whose pages display the non-triangulation survivors as
  filled n-gons with the chords drawn in red, two per row.
- Auto-opens the PDF in the default viewer.

The PDF is saved to:
    computations/visualizing nonlocals that refuse to die/
        nonlocal n=N (step1 survivors).pdf
"""

import math
import os
import subprocess
import sys
import time
from itertools import combinations_with_replacement

# --- plotting backend.  "Agg" means "no live window" -- we will write a PDF
# and open it in the OS PDF viewer instead, which is more portable than
# matplotlib's interactive windows. ----------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Polygon as MplPolygon


# ============================================================================
# 1.  COMBINATORICS OF CHORDS AND ZONES
# ============================================================================

def normalize(i, j, n):
    """
    Return the canonical name of the chord between vertices i and j of the
    n-gon, as a pair (a, b) with a < b and both in {1, ..., n}.

    We allow i and j to be any integers (negative, > n, etc.) and reduce
    them cyclically modulo n -- this lets the rest of the code write things
    like normalize(r + 1, r - 1, n) without worrying about wrap-around.
    """
    i = ((i - 1) % n) + 1     # bring i into {1, ..., n}
    j = ((j - 1) % n) + 1     # bring j into {1, ..., n}
    if i > j:
        i, j = j, i           # canonical order: smaller index first
    return (i, j)


def all_chords(n):
    """
    All planar chords of the convex n-gon -- i.e. all pairs (i, j) with
    1 <= i < j <= n that are NOT polygon edges.

    A pair is an edge iff j - i == 1 (consecutive vertices) or (i, j) ==
    (1, n) (the wrap-around edge from vertex n back to vertex 1).

    For an n-gon there are n(n-3)/2 chords.  Concretely:
        n = 4 -> 2 chords    (the two diagonals of the square)
        n = 5 -> 5 chords    (the pentagon's diagonals)
        n = 6 -> 9 chords
        n = 7 -> 14 chords, etc.
    """
    return [(i, j) for i in range(1, n + 1)
                   for j in range(i + 1, n + 1)
                   if j - i >= 2 and not (i == 1 and j == n)]


def zone_structure(r, n):
    """
    For zone Z_{r, r+2} return the data the kill criterion needs:

        special     -- the chord X_{r, r+2} sent to infinity in the limit
        bare        -- the chord X_{r+1, r-1} (cyclic) which equals -X_{r,r+2};
                       its 1/X factor also drops at leading order
        pairs       -- a list of (companion, substitute) chord-pairs
                       (X_{r, k}, X_{r+1, k}) for k = r+3, ..., r+n-2

    Each pair comes from one of the n-4 substitutions
        X_{r+1, k} = X_{r, k} - X_{r, r+2}
    that the 1-zero condition c_{r, k-1} = 0 enforces.  In the limit we
    must hold AT MOST ONE side of each pair finite (the other side goes
    to infinity), so a monomial whose multiset contains both elements of
    some pair has no consistent surviving limit -- it cannot be killed
    via this zone's Step-1 argument.
    """
    special = normalize(r,     r + 2, n)
    bare    = normalize(r + 1, r - 1, n)        # cyclic; non-adjacent chord

    pairs = []
    for offset in range(3, n - 1):              # offsets 3, 4, ..., n-2
        k = ((r - 1 + offset) % n) + 1          # cyclic position
        companion  = normalize(r,     k, n)     # X_{r, k}      (free side)
        substitute = normalize(r + 1, k, n)     # X_{r+1, k}    (substituted side)
        pairs.append((companion, substitute))

    return special, bare, pairs


def killable_at_zone(ms_set, structure):
    """
    Decide whether the multiset (passed in as a Python set of chord names
    for fast membership lookup) is killed by Step-1 at the zone described
    by `structure = (special, bare, pairs)`.

    The two conditions reproduced from the proof:

        (A)  multiset must NOT contain the special or the bare;
        (B)  for every (companion, substitute) pair, multiset must NOT
             contain BOTH chords of the pair.

    If both conditions hold, we can pick a sub-limit of the Z_{r,r+2}
    kinematics (one binary choice per pair: hold companion arbitrary, OR
    hold substitute arbitrary) consistent with the multiset, and the
    leading-order Laurent argument kills its coefficient.
    """
    special, bare, pairs = structure
    if special in ms_set or bare in ms_set:     # condition (A) fails
        return False
    for companion, substitute in pairs:         # condition (B) check
        if companion in ms_set and substitute in ms_set:
            return False
    return True


# ============================================================================
# 2.  GEOMETRY: CHORD CROSSING AND TRIANGULATIONS
# ============================================================================

def chords_cross(c1, c2):
    """
    Two chords (a, b) and (c, d) of the convex n-gon (with vertices in the
    cyclic order 1, 2, ..., n) cross in the interior iff their four
    endpoints alternate around the boundary -- equivalently, exactly one
    endpoint of (c, d) lies STRICTLY BETWEEN a and b (in the linear order
    a < ... < b), the other lies strictly outside.

    Sharing an endpoint => no crossing (chords meet at a vertex).
    """
    a, b = c1
    c, d = c2
    if {a, b} & {c, d}:                         # share a vertex
        return False
    inside = lambda v: a < v < b
    return inside(c) != inside(d)               # XOR


def is_triangulation(ms):
    """
    A multiset of chords is a triangulation of the n-gon iff
        (i)  all its chords are distinct (no double poles);
        (ii) no two chords cross.
    Then the chords carve the polygon into n-2 triangles -- this is what
    we call a "tree-amplitude diagram" (Catalan-many of them at each n).
    """
    if len(set(ms)) != len(ms):                 # repeated chord => double pole
        return False
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            if chords_cross(ms[i], ms[j]):
                return False
    return True


# ============================================================================
# 3.  DIAGRAMS: DRAW EACH SURVIVING MULTISET AS AN n-GON WITH RED CHORDS
# ============================================================================

def vertex_xy(i, n, R=1.0):
    """
    Cartesian coordinates of vertex i of a regular n-gon inscribed in a
    circle of radius R.  Vertex 1 sits at the top; vertices proceed
    clockwise (1, 2, 3, ...) -- matching the cyclic order used everywhere
    above.
    """
    angle = math.pi / 2 - 2 * math.pi * (i - 1) / n
    return (R * math.cos(angle), R * math.sin(angle))


def draw_diagram(ax, multiset, n):
    """
    Draw a single survivor on the supplied matplotlib axis:
        - filled black n-gon
        - each chord in the multiset rendered as a red line segment
        - if a chord appears with multiplicity > 1 (a double pole), draw
          parallel offset copies so the multiplicity is visible
        - print the coefficient label  a_{ {i1,j1}, {i2,j2}, ... }  on the
          right of the polygon, in the same notation used in your existing
          n=8 PDF.
    """
    coords = [vertex_xy(i, n) for i in range(1, n + 1)]
    ax.add_patch(MplPolygon(coords, closed=True,
                            facecolor='black', edgecolor='black'))

    # count how many times each distinct chord appears (handle double poles)
    chord_counts = {}
    for chord in multiset:
        chord_counts[chord] = chord_counts.get(chord, 0) + 1

    for chord, mult in chord_counts.items():
        i, j = chord
        x1, y1 = vertex_xy(i, n)
        x2, y2 = vertex_xy(j, n)
        # for multiplicity m, draw m parallel copies offset perpendicular
        # to the chord direction so they are visually distinguishable
        for k in range(mult):
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy) or 1.0
            nx, ny = -dy / length, dx / length          # unit normal
            shift = 0.025 * (k - (mult - 1) / 2.0)      # symmetric offsets
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
    """
    Pack the survivor diagrams into multiple letter-size pages, two columns
    of four rows each (8 diagrams per page).  Returns a list of figures
    that the caller will dump into a PDF.
    """
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
        # blank out unused subplots
        axes_flat = axes.flatten() if rows * cols > 1 else [axes]
        for ax in axes_flat:
            ax.axis('off')
        for ax, ms in zip(axes_flat, page):
            draw_diagram(ax, ms, n)
        fig.tight_layout(rect=[0, 0, 0.95, 0.96])
        figs.append(fig)
    return figs


# ============================================================================
# 4.  MAIN: GET n, ENUMERATE, CLASSIFY, RENDER, OPEN
# ============================================================================

def parse_n():
    """
    Get n from the command line if provided, otherwise prompt interactively.
    Recommended cap is n <= 9 in pure Python; see WHEN-IT-BREAKS notes
    in the project README.
    """
    if len(sys.argv) > 1:
        try:
            return int(sys.argv[1])
        except ValueError:
            pass
    raw = input("Enter n (>= 4, recommended <= 9): ").strip()
    return int(raw)


def main():
    n = parse_n()
    if n < 4:
        sys.exit(f"n = {n} too small; need n >= 4.")
    N = n - 3                                   # multiset size = chords per monomial

    chords     = all_chords(n)
    structures = [zone_structure(r, n) for r in range(1, n + 1)]

    # --- print zone tables so the reader can match what the program is
    #     doing against the by-hand calculation in the notes ------------
    print(f"\nn = {n}, multiset size N = {N}, |chords| = {len(chords)}")
    print("Zone structures (special, bare, [companion <-> substituted pairs]):")
    for r, (sp, br, pr) in enumerate(structures, start=1):
        rp2 = normalize(r, r + 2, n)
        print(f"  Z_{rp2}:  sp={sp}  bare={br}  pairs={pr}")

    # --- enumerate every size-N multiset of chords and apply the kill
    #     criterion at every zone.  A multiset is "uncaught" iff NO zone
    #     can kill it. ----------------------------------------------------
    print(f"\nEnumerating multisets...")
    t0 = time.time()
    total = 0
    uncaught = []
    for ms in combinations_with_replacement(chords, N):
        total += 1
        ms_set = set(ms)                       # turn into a fast lookup table
        if not any(killable_at_zone(ms_set, s) for s in structures):
            uncaught.append(ms)
    dt = time.time() - t0

    # --- separate the survivors into triangulations (expected -- they are
    #     the tree amplitude itself) and non-triangulations (the coefficients
    #     that Step-1 cannot eliminate -- the "diagrams that refuse to die").
    tri    = [ms for ms in uncaught if     is_triangulation(ms)]
    nontri = [ms for ms in uncaught if not is_triangulation(ms)]

    print(f"  ({dt:.2f} s)")
    print(f"\nTotal multisets of size {N} from {len(chords)} chords: {total}")
    print(f"Uncaught by step-1 at any of {n} zones:                 {len(uncaught)}")
    print(f"  -- triangulations among uncaught:                       {len(tri)}")
    print(f"  -- non-triangulation uncaught (need higher layers):     {len(nontri)}")

    if nontri and len(nontri) <= 30:
        print(f"\nNon-triangulation uncaught multisets ({len(nontri)}):")
        for ms in nontri:
            print(f"  {ms}")
    elif nontri:
        print(f"\n(omitting print of {len(nontri)} multisets; see PDF)")

    # --- render to PDF in the project's diagram folder -----------------
    here     = os.path.dirname(os.path.abspath(__file__))
    out_dir  = os.path.join(here, "outputs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"nonlocal n={n} (step1 survivors).pdf")

    print(f"\nRendering {len(nontri)} diagram(s)...")
    figs = build_figures(nontri, n,
        title=f"n = {n}: non-triangulation multisets uncaught by Step-1")
    with PdfPages(out_path) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Wrote: {out_path}")

    # --- auto-open the PDF in the OS default viewer ---------------------
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


if __name__ == "__main__":
    main()
