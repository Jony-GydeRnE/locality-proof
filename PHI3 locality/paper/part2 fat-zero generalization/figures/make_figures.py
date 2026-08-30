#!/usr/bin/env python3
"""
Figure generator for "Part II: the fat-zero generalisation".
Produces clean vector (PDF) n-gon diagrams: the k-zero geometry, the
(k+2)-gon-over-(k+1)-gon survival, the held-finite supports of regimes
A and B, the per-column slots of regime C, and two survivor pictures.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

OUT = os.path.dirname(os.path.abspath(__file__))

# ----- colours -----
C_EDGE   = "0.12"
C_NEAR   = "#f2c94c"   # near-block shade
C_FAR    = "#9bc4e2"   # far-block shade
C_ENCL   = "#1f9d55"   # enclosing chord (length-k)
C_SPECIAL= "#d64541"   # special chord
C_BARE   = "#2e5fa3"   # bare chord
C_FREE   = "#e8820e"   # held-finite (free) bridges / survivor chords
C_CROSS  = "#9b51e0"   # crossing chord (P1)


def vpos(n):
    ang = np.pi / 2 - 2 * np.pi * np.arange(n) / n   # vertex i (0-indexed), clockwise from top
    return np.cos(ang), np.sin(ang)


def setup(ax, pad=1.34):
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-pad, pad)
    ax.set_ylim(-pad, pad)


def ngon(ax, n, labels=True):
    x, y = vpos(n)
    for i in range(n):
        j = (i + 1) % n
        ax.plot([x[i], x[j]], [y[i], y[j]], color=C_EDGE, lw=2.0, zorder=2)
    ax.scatter(x, y, s=26, color=C_EDGE, zorder=5)
    if labels:
        for i in range(n):
            ax.text(1.17 * x[i], 1.17 * y[i], str(i + 1),
                    ha="center", va="center", fontsize=10.5, zorder=6)


def chord(ax, n, a, b, **kw):
    x, y = vpos(n)
    ax.plot([x[a - 1], x[b - 1]], [y[a - 1], y[b - 1]], zorder=4, **kw)


def shade(ax, n, verts, color, alpha=0.5):
    x, y = vpos(n)
    pts = [(x[v - 1], y[v - 1]) for v in verts]
    ax.add_patch(MplPolygon(pts, closed=True, facecolor=color,
                            edgecolor="none", alpha=alpha, zorder=1))


def panel_tag(ax, s):
    ax.text(0.5, -0.06, s, transform=ax.transAxes, ha="center", va="top",
            fontsize=10)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


# ---------------------------------------------------------------- Fig 1
def fig_geometry():
    n, k = 10, 3
    fig, ax = plt.subplots(figsize=(3.4, 3.4))
    setup(ax)
    near = list(range(1, k + 2))            # 1..k+1
    far = list(range(k + 2, n + 1))         # k+2..n
    shade(ax, n, near, C_NEAR, 0.55)
    shade(ax, n, far, C_FAR, 0.30)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=3.0)                     # enclosing X_{1,k+1}
    chord(ax, n, 1, k + 2, color=C_SPECIAL, lw=2.2, ls=(0, (5, 2)))  # special X_{1,k+2}
    chord(ax, n, k + 1, n, color=C_BARE, lw=2.2, ls=(0, (5, 2)))     # bare X_{k+1,n}
    save(fig, "geometry.pdf")


# ---------------------------------------------------------------- Fig 2
def fig_kplus2gon():
    n, k = 10, 3
    m = 8  # interior far vertex for the column-pair case
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.4))
    for ax in axs:
        setup(ax)
    # (a) close on the special side: (k+2)-gon {1,..,k+2}
    ax = axs[0]
    shade(ax, n, list(range(1, k + 3)), C_NEAR, 0.5)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=3.0)
    chord(ax, n, 1, k + 2, color=C_SPECIAL, lw=2.4)
    panel_tag(ax, r"(a)  via special $X_{1,k+2}$")
    # (b) close on the bare side: (k+2)-gon {n,1,..,k+1}
    ax = axs[1]
    shade(ax, n, [n] + list(range(1, k + 2)), C_NEAR, 0.5)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=3.0)
    chord(ax, n, k + 1, n, color=C_BARE, lw=2.4)
    panel_tag(ax, r"(b)  via bare $X_{k+1,n}$")
    # (c) column-pair: (k+2)-gon {1,..,k+1, m} closed by {X_{1,m}, X_{k+1,m}}
    ax = axs[2]
    shade(ax, n, list(range(1, k + 2)) + [m], C_NEAR, 0.5)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=3.0)
    chord(ax, n, 1, m, color=C_FREE, lw=2.4)
    chord(ax, n, k + 1, m, color=C_FREE, lw=2.4)
    panel_tag(ax, r"(c)  via column pair $\{X_{1,m},\,X_{k+1,m}\}$")
    save(fig, "kplus2gon.pdf")


# ---------------------------------------------------------------- Fig 3 (NEW)
# Recursive stacking: as k grows, the (k+2)-gon at the k-zero is itself the
# (k+1)-gon at the (k+1)-zero -- so the polygons build on top of each other.
def fig_stacking():
    n = 10
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.4))
    for ax, k in zip(axs, [1, 2, 3]):
        setup(ax)
        # the (k+1)-gon over k edges = vertices 1..k+1
        shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.45)
        # closing one more vertex -> (k+2)-gon  = vertices 1..k+2
        shade(ax, n, list(range(1, k + 3)), C_NEAR, 0.20)
        ngon(ax, n)
        # enclosing chord of the (k+1)-gon
        if k >= 2:
            chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.6)
        # the new chord that closes into the (k+2)-gon (here: via the special)
        chord(ax, n, 1, k + 2, color=C_SPECIAL, lw=2.6)
        panel_tag(ax, rf"$k = {k}$:  ({k+1})-gon $\Rightarrow$ ({k+2})-gon")
    save(fig, "stacking.pdf")


# ---------------------------------------------------------------- Fig: survivors
# Three geometric survival modes: P2 (build a (k+2)-gon via special / bare /
# column pair) and P1 (a chord crosses the enclosing chord).
def fig_survivors():
    n, k = 8, 2  # k=2 zero based at {1,2,3}, enclosing X_{1,3}
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.4))
    # (a) P2 via special: (k+2)-gon {1,2,3,4} closed by special X_{1,4}
    ax = axs[0]
    setup(ax)
    shade(ax, n, list(range(1, k + 3)), C_NEAR, 0.4)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.2)
    chord(ax, n, 1, k + 2, color=C_SPECIAL, lw=2.6)
    chord(ax, n, 5, 7, color=C_FREE, lw=2.0)
    chord(ax, n, 5, 8, color=C_FREE, lw=2.0)
    panel_tag(ax, r"P2 via special $X_{1,4}$")
    # (b) P2 via column pair {X_{1,6}, X_{3,6}}: (k+2)-gon {1,2,3,6}
    ax = axs[1]
    setup(ax)
    shade(ax, n, list(range(1, k + 2)) + [6], C_NEAR, 0.4)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.2)
    chord(ax, n, 1, 6, color=C_FREE, lw=2.6)
    chord(ax, n, k + 1, 6, color=C_FREE, lw=2.6)
    chord(ax, n, 5, 7, color=C_FREE, lw=2.0)
    panel_tag(ax, r"P2 via column pair $\{X_{1,6},\,X_{3,6}\}$")
    # (c) P1: a chord crosses the enclosing chord X_{1,3}
    ax = axs[2]
    setup(ax)
    shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.4)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.2)
    chord(ax, n, 2, 6, color=C_CROSS, lw=2.6)          # crosses X_{1,3}
    chord(ax, n, 5, 7, color=C_FREE, lw=2.0)
    chord(ax, n, 4, 8, color=C_FREE, lw=2.0)
    panel_tag(ax, r"P1: $X_{2,6}$ crosses enclosing $X_{1,3}$")
    save(fig, "survivors.pdf")


# ---------------------------------------------------------------- Fig: per-column choice
# At a k-zero, each far vertex l carries a "column" of k+1 chords coupled by
# the chord identities. The asymptotic analysis shows that under the kill limit
# we can hold finite exactly one of three slots per column. (Notes 24-25.)
def fig_column_rule():
    n, k, ell = 8, 2, 6
    fig, axs = plt.subplots(1, 3, figsize=(8.0, 2.8))
    for ax, (chords, tag) in zip(axs, [
        ([(1, ell)],                            r"hold $X_{1,\ell}$ finite"),
        ([(k + 1, ell)],                        r"hold $X_{k+1,\ell}$ finite"),
        ([(i, ell) for i in range(2, k + 1)],   r"hold $\{X_{2,\ell},\dots,X_{k,\ell}\}$ finite"),
    ]):
        setup(ax)
        shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.4)
        ngon(ax, n)
        chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.2)
        chord(ax, n, 1, k + 2, color=C_SPECIAL, lw=1.2, ls=(0, (3, 2)), alpha=0.5)
        chord(ax, n, k + 1, n, color=C_BARE, lw=1.2, ls=(0, (3, 2)), alpha=0.5)
        for (a, b) in chords:
            chord(ax, n, a, b, color=C_FREE, lw=2.4)
        panel_tag(ax, tag)
    save(fig, "column_rule.pdf")


# ---------------------------------------------------------------- Fig: surviving multisets gallery
# Concrete surviving multisets at n=7, k=2 (the 2-zero based at vertices {1,2,3}).
# Each shows a multiset of size n-3 = 4 that contains an incompatible chord pair,
# so the rate equations from the chord identities have no consistent solution.
def fig_survivor_gallery():
    n, k = 7, 2
    fig, axs = plt.subplots(1, 4, figsize=(11.0, 2.8))

    def draw_one(ax, M_pairs, tag, highlight_pair=None):
        setup(ax)
        shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.35)
        ngon(ax, n)
        chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.0)
        for (a, b) in M_pairs:
            col = C_FREE
            lw = 2.0
            if highlight_pair and (a, b) in highlight_pair:
                col = C_CROSS
                lw = 2.6
            chord(ax, n, a, b, color=col, lw=lw)
        panel_tag(ax, tag)

    # (a) special X_{1,4} present (always survives)
    draw_one(axs[0],
             [(1, 4), (1, 3), (4, 6), (5, 7)],
             r"(a) special $X_{1,4}\in M$")
    # (b) bare X_{3,7} present
    draw_one(axs[1],
             [(3, 7), (1, 3), (4, 6), (5, 7)],
             r"(b) bare $X_{3,7}\in M$")
    # (c) column pair {X_{1,5}, X_{3,5}} (force alpha_5=0 AND alpha_5=1, impossible)
    draw_one(axs[2],
             [(1, 5), (3, 5), (1, 3), (4, 6)],
             r"(c) column pair $\{X_{1,5},X_{3,5}\}$",
             highlight_pair={(1, 5), (3, 5)})
    # (d) offset pair {X_{2,4}, X_{2,7}} (force beta_2=0 AND beta_2=1, impossible)
    draw_one(axs[3],
             [(2, 4), (2, 7), (1, 3), (4, 6)],
             r"(d) offset pair $\{X_{2,4},X_{2,7}\}$",
             highlight_pair={(2, 4), (2, 7)})
    save(fig, "survivor_gallery.pdf")


# ---------------------------------------------------------------- Fig: fan of enclosing chords
# From a single base vertex i, the enclosing chords X_{i, i+k} for k = 1..floor(n/2)
# form a fan. Surviving every k-zero based at vertex i via P1 alone requires a
# crossing chord for each enclosing chord -- ~n/2 crossings in total.
def fig_fan():
    n = 10
    fig, ax = plt.subplots(figsize=(3.6, 3.6))
    setup(ax)
    ngon(ax, n)
    # Draw the fan from vertex 1: X_{1, 1+k} for k = 2, 3, ..., floor(n/2)
    # (k=1 gives the leg X_{1,2}, k=floor(n/2) gives the longest chord)
    fan_colors = ["#1f9d55", "#1f9d55", "#1f9d55", "#1f9d55"]
    for k, c in zip(range(2, n // 2 + 1), fan_colors):
        chord(ax, n, 1, 1 + k, color=c, lw=1.8)
    # Add example crossings (purple): chords that cross each enclosing chord
    # These would be the simultaneous P1 crossings needed.
    chord(ax, n, 2, 10, color=C_CROSS, lw=1.6, ls=(0, (3, 2)))  # crosses X_{1,3}
    chord(ax, n, 3, 9, color=C_CROSS, lw=1.6, ls=(0, (3, 2)))   # crosses X_{1,4}
    chord(ax, n, 4, 8, color=C_CROSS, lw=1.6, ls=(0, (3, 2)))   # crosses X_{1,5}
    save(fig, "fan.pdf")


if __name__ == "__main__":
    fig_geometry()
    fig_kplus2gon()
    fig_stacking()
    fig_survivors()
    fig_column_rule()
    fig_survivor_gallery()
    fig_fan()
    print("figures written to", OUT)
