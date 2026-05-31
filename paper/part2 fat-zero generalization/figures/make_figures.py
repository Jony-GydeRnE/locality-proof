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


# ---------------------------------------------------------------- Fig 3 & 4
def _regime(ax, n, k, bridges, dashed_special_bare=True):
    setup(ax)
    shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.45)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.6)
    if dashed_special_bare:
        chord(ax, n, 1, k + 2, color=C_SPECIAL, lw=1.6, ls=(0, (3, 2)), alpha=0.7)
        chord(ax, n, k + 1, n, color=C_BARE, lw=1.6, ls=(0, (3, 2)), alpha=0.7)
    for (a, b) in bridges:
        chord(ax, n, a, b, color=C_FREE, lw=2.2)


def fig_regimeA():
    n, k = 9, 2
    # A : {X_{i,k+2}: i=2..k} U {X_{1,j}: j=k+3..n-1}
    bridges = [(i, k + 2) for i in range(2, k + 1)] + [(1, j) for j in range(k + 3, n)]
    fig, ax = plt.subplots(figsize=(3.4, 3.4))
    _regime(ax, n, k, bridges)
    save(fig, "regimeA.pdf")


def fig_regimeB():
    n, k = 9, 2
    # B : {X_{i,n}: i=2..k} U {X_{k+1,j}: j=k+3..n-1}
    bridges = [(i, n) for i in range(2, k + 1)] + [(k + 1, j) for j in range(k + 3, n)]
    fig, ax = plt.subplots(figsize=(3.4, 3.4))
    _regime(ax, n, k, bridges)
    save(fig, "regimeB.pdf")


# ---------------------------------------------------------------- Fig 5
def fig_criterionC():
    n, k, l = 8, 2, 6
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
    slots = [
        ([(1, l)],                              r"slot $X_{1,\ell}$"),
        ([(k + 1, l)],                          r"slot $X_{k+1,\ell}$"),
        ([(i, l) for i in range(2, k + 1)],     r"slot $\{X_{2,\ell},\dots,X_{k,\ell}\}$"),
    ]
    for ax, (bridges, tag) in zip(axs, slots):
        setup(ax)
        shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.4)
        ngon(ax, n)
        chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.2)
        chord(ax, n, 1, k + 2, color=C_SPECIAL, lw=1.3, ls=(0, (3, 2)), alpha=0.55)
        chord(ax, n, k + 1, n, color=C_BARE, lw=1.3, ls=(0, (3, 2)), alpha=0.55)
        for (a, b) in bridges:
            chord(ax, n, a, b, color=C_FREE, lw=2.4)
        panel_tag(ax, tag)
    save(fig, "criterionC.pdf")


# ---------------------------------------------------------------- Fig 6
def fig_survivors():
    n, k = 8, 2  # bigger so a column-pair example fits
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.4))
    # (a) special present -> (k+2)-gon completion {1,..,k+2}
    ax = axs[0]
    setup(ax)
    shade(ax, n, list(range(1, k + 3)), C_NEAR, 0.4)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.2)
    chord(ax, n, 1, k + 2, color=C_SPECIAL, lw=2.6)
    chord(ax, n, 5, 7, color=C_FREE, lw=2.0)
    chord(ax, n, 5, 8, color=C_FREE, lw=2.0)
    panel_tag(ax, r"(a) special $X_{1,4}$ present")
    # (b) column-pair survival: X_{1,m} AND X_{k+1,m} for m=6 -> rate incompat
    ax = axs[1]
    setup(ax)
    shade(ax, n, list(range(1, k + 2)) + [6], C_NEAR, 0.4)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.2)
    chord(ax, n, 1, 6, color=C_FREE, lw=2.6)
    chord(ax, n, k + 1, 6, color=C_FREE, lw=2.6)
    chord(ax, n, 5, 7, color=C_FREE, lw=2.0)
    panel_tag(ax, r"(b) column pair $\{X_{1,6},X_{3,6}\}$")
    # (c) column-(k+2) x column-n pair survival: X_{i,k+2} AND X_{i,n}, same i
    ax = axs[2]
    setup(ax)
    shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.4)
    ngon(ax, n)
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.2)
    chord(ax, n, 2, k + 2, color=C_FREE, lw=2.6)      # X_{2, k+2}
    chord(ax, n, 2, n, color=C_FREE, lw=2.6)          # X_{2, n}
    chord(ax, n, 5, 7, color=C_FREE, lw=2.0)
    panel_tag(ax, r"(c) offset pair $\{X_{2,4},X_{2,8}\}$")
    save(fig, "survivors.pdf")


if __name__ == "__main__":
    fig_geometry()
    fig_kplus2gon()
    fig_regimeA()
    fig_regimeB()
    fig_criterionC()
    fig_survivors()
    print("figures written to", OUT)
