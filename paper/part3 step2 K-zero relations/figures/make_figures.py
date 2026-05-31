#!/usr/bin/env python3
"""
Figures for Part III: Step-2 K-zero relations.
Diagrammatic identities (multiset coefficients drawn on the n-gon) and the
boundary-propagator structure.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

OUT = os.path.dirname(os.path.abspath(__file__))

C_EDGE   = "0.12"
C_NEAR   = "#f2c94c"
C_FAR    = "#9bc4e2"
C_ENCL   = "#1f9d55"
C_SPECIAL= "#d64541"
C_BARE   = "#2e5fa3"
C_FREE   = "#e8820e"
C_BP     = "#9b51e0"     # boundary-propagator highlight
C_BG     = "#cfe8c9"     # green "blob" for the rest of the diagram


def vpos(n):
    ang = np.pi / 2 - 2 * np.pi * np.arange(n) / n
    return np.cos(ang), np.sin(ang)


def setup(ax, pad=1.34):
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-pad, pad); ax.set_ylim(-pad, pad)


def ngon(ax, n, labels=True):
    x, y = vpos(n)
    for i in range(n):
        j = (i + 1) % n
        ax.plot([x[i], x[j]], [y[i], y[j]], color=C_EDGE, lw=2.0, zorder=2)
    ax.scatter(x, y, s=24, color=C_EDGE, zorder=5)
    if labels:
        for i in range(n):
            ax.text(1.17 * x[i], 1.17 * y[i], str(i + 1),
                    ha="center", va="center", fontsize=10, zorder=6)


def chord(ax, n, a, b, **kw):
    x, y = vpos(n)
    ax.plot([x[a - 1], x[b - 1]], [y[a - 1], y[b - 1]], zorder=4, **kw)


def shade(ax, n, verts, color, alpha=0.5):
    x, y = vpos(n)
    pts = [(x[v - 1], y[v - 1]) for v in verts]
    ax.add_patch(MplPolygon(pts, closed=True, facecolor=color,
                            edgecolor="none", alpha=alpha, zorder=1))


def panel_tag(ax, s, dy=-0.06):
    ax.text(0.5, dy, s, transform=ax.transAxes, ha="center", va="top", fontsize=10)


def op_text(ax, s, x, y, **kw):
    ax.text(x, y, s, ha="center", va="center", fontsize=18, **kw)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


# ------------------------------------------------------------- Fig 1
def fig_bp_structure():
    """A schematic 2-zero (n=8): special, bare, born-free 'blob', 2-BP example."""
    n, k = 8, 2
    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    setup(ax)
    # near block shaded
    shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.5)
    # blob over far block
    shade(ax, n, list(range(k + 2, n + 1)), C_BG, 0.7)
    ngon(ax, n)
    # enclosing
    chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.4)
    # special X_{1,4} (highlighted as 'BP')
    chord(ax, n, 1, k + 2, color=C_SPECIAL, lw=2.4)
    # bare X_{3,n}
    chord(ax, n, k + 1, n, color=C_BARE, lw=2.4)
    # one BP example bridge: X_{1,6}
    chord(ax, n, 1, 6, color=C_BP, lw=2.0, ls=(0, (4, 2)))
    save(fig, "bp_structure.pdf")


# ------------------------------------------------------------- Fig 2
def fig_binomial_trace():
    """The 2-BP 'binomial-trace' identity for the 2-zero (k=2, n=8):
       a_{14,14} - a_{14,3n} + a_{3n,3n} = 0
       drawn as three diagrams summing to 0 with signs.
    """
    n, k = 8, 2
    fig, axs = plt.subplots(1, 5, figsize=(11.6, 2.6),
                            gridspec_kw=dict(width_ratios=[1, 0.18, 1, 0.18, 1]))
    # subplot positions: 0 diagram, 1 sign, 2 diagram, 3 sign, 4 diagram
    diagram_axes = [axs[0], axs[2], axs[4]]
    sign_axes    = [axs[1], axs[3]]

    def base_panel(ax, chords, tag):
        setup(ax)
        shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.35)
        shade(ax, n, list(range(k + 2, n + 1)), C_BG, 0.6)
        ngon(ax, n)
        chord(ax, n, 1, k + 1, color=C_ENCL, lw=2.0)
        for (a, b, col) in chords:
            chord(ax, n, a, b, color=col, lw=2.6)
        panel_tag(ax, tag)

    # (a)  a_{14,14}
    base_panel(diagram_axes[0],
               [(1, k + 2, C_SPECIAL), (1, k + 2, C_SPECIAL)],
               r"$a_{14,14}$")
    # minus
    sign_axes[0].axis("off"); op_text(sign_axes[0], "$-$", 0.5, 0.5)
    # (b)  a_{14,3n}
    base_panel(diagram_axes[1],
               [(1, k + 2, C_SPECIAL), (k + 1, n, C_BARE)],
               r"$a_{14,3n}$")
    sign_axes[1].axis("off"); op_text(sign_axes[1], "$+$", 0.5, 0.5)
    # (c)  a_{3n,3n}
    base_panel(diagram_axes[2],
               [(k + 1, n, C_BARE), (k + 1, n, C_BARE)],
               r"$a_{3n,3n}\ =\ 0$")
    save(fig, "binomial_trace.pdf")


# ------------------------------------------------------------- Fig 3
def fig_column_relation():
    """The 2-zero 2-BP column relation at j=5 (k=2, n=8):
       a_{14,15} + a_{14,25} + a_{14,35} - a_{3n,15} - a_{3n,25} - a_{3n,35} = 0
       (Six diagrams: special paired with column-5 chords, minus bare paired
        with column-5 chords.)
    """
    n, k, j = 8, 2, 5
    fig, axs = plt.subplots(1, 11, figsize=(13.6, 2.0),
                            gridspec_kw=dict(width_ratios=[1, 0.15] * 5 + [1]))

    def base_panel(ax, chords, tag):
        setup(ax)
        shade(ax, n, list(range(1, k + 2)), C_NEAR, 0.3)
        shade(ax, n, list(range(k + 2, n + 1)), C_BG, 0.55)
        ngon(ax, n)
        chord(ax, n, 1, k + 1, color=C_ENCL, lw=1.6)
        for (a, b, col) in chords:
            chord(ax, n, a, b, color=col, lw=2.2)
        panel_tag(ax, tag, dy=-0.10)

    s_chord = (1, k + 2, C_SPECIAL)
    b_chord = (k + 1, n, C_BARE)
    col_chords = [(1, j), (2, j), (k + 1, j)]
    labels = [r"$a_{14,15}$", r"$a_{14,25}$", r"$a_{14,35}$",
              r"$a_{3n,15}$", r"$a_{3n,25}$", r"$a_{3n,35}$"]
    # pluses then minuses
    signs = ["$+$", "$+$", "$-$", "$-$", "$-$"]

    # panels 0,2,4,6,8,10 ; sign axes 1,3,5,7,9
    diagram_axes = [axs[i] for i in [0, 2, 4, 6, 8, 10]]
    sign_axes    = [axs[i] for i in [1, 3, 5, 7, 9]]

    for i in range(3):
        c = (col_chords[i][0], col_chords[i][1], C_FREE)
        base_panel(diagram_axes[i], [s_chord, c], labels[i])
    for i in range(3):
        c = (col_chords[i][0], col_chords[i][1], C_FREE)
        base_panel(diagram_axes[3 + i], [b_chord, c], labels[3 + i])
    for k_, s_ax in enumerate(sign_axes):
        s_ax.axis("off"); op_text(s_ax, signs[k_], 0.5, 0.5)
    diagram_axes[-1].text(1.20, 0.10, "$=\\ 0$",
                          transform=diagram_axes[-1].transAxes,
                          ha="left", va="center", fontsize=14)
    save(fig, "column_relation.pdf")


if __name__ == "__main__":
    fig_bp_structure()
    fig_binomial_trace()
    fig_column_relation()
    print("figures written to", OUT)
