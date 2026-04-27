"""
Generate the (a)+(b) worldline figure for §2.4 Part 1 of geometric_story.tex.
Output: figures/worldline_figure.pdf  (referenced via \\includegraphics)

Key design points:
  * Common starting event START = (x0, t0) = (1.0, 1.5) -- shifted off the
    axes so the "at rest" line is clearly visible (not hugging the t-axis)
    and so the curves can asymptotically approach the lightcone from above.
  * The lightcone is the lightcone OF THE STARTING EVENT, not of the origin:
    its right branch is t = t0 + (x - x0) for x >= x0, drawn dashed.
  * Worldlines are correct relativistic hyperbolas with constant proper
    acceleration a, parameterized as
        x(tau) = x0 + (1/a)(cosh(a*tau) - 1)
        t(tau) = t0 + (1/a) sinh(a*tau)
    so they pass through (x0, t0) at rest (dx/dt|_{tau=0} = 0) and asymptote
    to slope dt/dx = 1 (i.e. v = c) from above. They stay strictly INSIDE
    the lightcone for all tau > 0.
  * Slow particle (small a): asymptote is well above the lightcone -> curve
    stays close to the t-axis (more vertical worldline).
  * Fast particle (large a): asymptote is just above the lightcone -> curve
    tilts sharply right, hugging the lightcone.
  * Panel (b) shares the same starting event; the closed-loop arc rises into
    the future, peaks, and descends back to a "return" point at smaller t,
    illustrating that returning to the past requires v > c.

Run: python3 figures/make_worldline_figure.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "mathtext.fontset": "cm",
})

fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0))

# common starting event for both panels
x0, t0 = 1.0, 1.5
xmax, ymax = 5.0, 4.7

WORLDLINE = "#1f3a93"
LOOP      = "#a02020"

# ---------------------------------------------------------- panel (a) ----
ax = axes[0]

# axes
ax.annotate("", xy=(xmax, 0), xytext=(-0.3, 0),
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
ax.annotate("", xy=(0, ymax), xytext=(0, -0.3),
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
ax.text(xmax + 0.05, -0.05, "space", ha="left", va="top")
ax.text(-0.1, ymax + 0.05, "time", ha="right", va="bottom")

# lightcone of the starting event (right branch only — the relevant one)
ax.plot([x0, xmax - 0.2], [t0, t0 + (xmax - 0.2 - x0)],
        "--", color="gray", lw=0.9)
ax.text(xmax - 0.05, t0 + (xmax - 0.2 - x0),
        r"$v=c$", color="gray", fontsize=9, ha="left", va="bottom")

# starting event marker
ax.plot([x0], [t0], "o", color="black", markersize=4)
ax.text(x0 - 0.1, t0 - 0.05, r"start", color="black", fontsize=9,
        ha="right", va="top")

# (1) "at rest" — vertical worldline rising from start
ax.plot([x0, x0], [t0, ymax - 0.25], color=WORLDLINE, lw=1.7)
ax.text(x0 - 0.12, (t0 + ymax) / 2, "at rest",
        color=WORLDLINE, fontsize=9, ha="right", va="center")

# (2) slow particle — small acceleration, asymptote far above lightcone
a_slow = 0.45
tau = np.linspace(0, 6, 400)
x_slow = x0 + (1.0 / a_slow) * (np.cosh(a_slow * tau) - 1)
t_slow = t0 + (1.0 / a_slow) * np.sinh(a_slow * tau)
mask = (x_slow <= xmax - 0.3) & (t_slow <= ymax - 0.25)
ax.plot(x_slow[mask], t_slow[mask], color=WORLDLINE, lw=1.7)
i = int(0.85 * mask.sum())
ax.text(x_slow[mask][i] - 0.15, t_slow[mask][i] + 0.02, "slow",
        color=WORLDLINE, fontsize=9, ha="right", va="bottom")

# (3) fast particle — large acceleration, asymptote close to lightcone
a_fast = 1.6
tau = np.linspace(0, 3, 400)
x_fast = x0 + (1.0 / a_fast) * (np.cosh(a_fast * tau) - 1)
t_fast = t0 + (1.0 / a_fast) * np.sinh(a_fast * tau)
mask = (x_fast <= xmax - 0.3) & (t_fast <= ymax - 0.25)
ax.plot(x_fast[mask], t_fast[mask], color=WORLDLINE, lw=1.7)
i = int(0.7 * mask.sum())
ax.annotate(r"fast: $v \to c$",
            xy=(x_fast[mask][i], t_fast[mask][i]),
            xytext=(x_fast[mask][i] + 0.5, t_fast[mask][i] - 0.7),
            color=WORLDLINE, fontsize=9,
            arrowprops=dict(arrowstyle="-", color=WORLDLINE, lw=0.6))

ax.set_xlim(-0.6, xmax + 0.3)
ax.set_ylim(-0.6, ymax + 0.3)
ax.set_aspect("equal")
ax.axis("off")
ax.text((xmax) / 2, -0.95,
        "(a) physical worldlines stay inside the lightcone",
        ha="center", va="top", fontsize=9)

# ---------------------------------------------------------- panel (b) ----
ax = axes[1]

ax.annotate("", xy=(xmax, 0), xytext=(-0.3, 0),
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
ax.annotate("", xy=(0, ymax), xytext=(0, -0.3),
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
ax.text(xmax + 0.05, -0.05, "space", ha="left", va="top")
ax.text(-0.1, ymax + 0.05, "time", ha="right", va="bottom")

# lightcone of starting event (right branch)
ax.plot([x0, xmax - 0.2], [t0, t0 + (xmax - 0.2 - x0)],
        "--", color="gray", lw=0.9)

# starting event marker (same as panel a)
ax.plot([x0], [t0], "o", color=LOOP, markersize=4)
ax.text(x0 - 0.1, t0 - 0.05, "start", color=LOOP, fontsize=9,
        ha="right", va="top")

# closed-timelike-curve arc: rises from start into future, peaks, then
# descends back DOWN in t to a "return" point with t_ret < t0
t_start = t0
t_return = t0 - 0.6  # below the starting time -> in the past
x_return = x0 + 3.4
# parametric ellipse-like arc going from (x0,t_start) to (x_return,t_return)
theta = np.linspace(np.pi, 0, 200)
cx = (x0 + x_return) / 2.0
cy = (t_start + t_return) / 2.0
rx = (x_return - x0) / 2.0
ry = 1.9
x_arc = cx + rx * np.cos(theta)
y_arc = cy + ry * np.sin(theta)
ax.plot(x_arc, y_arc, color=LOOP, lw=1.9)
# arrowhead near the descending tail
ax.annotate("", xy=(x_arc[-2], y_arc[-2]), xytext=(x_arc[-15], y_arc[-15]),
            arrowprops=dict(arrowstyle="->", color=LOOP, lw=1.7))

# return point
ax.plot([x_return], [t_return], "o", color=LOOP, markersize=4)
ax.text(x_return + 0.1, t_return - 0.05, "return",
        color=LOOP, fontsize=9, ha="left", va="top")

# annotation of the descending portion
ax.text(x_return + 0.2, (t_start + t_return) / 2 + 0.5,
        "curve goes\n" + r"$\mathit{down}$ in $t$" + "\n" + r"($v > c$)",
        color=LOOP, fontsize=9, ha="left", va="center")

ax.set_xlim(-0.6, xmax + 0.3)
ax.set_ylim(-0.6, ymax + 0.3)
ax.set_aspect("equal")
ax.axis("off")
ax.text((xmax) / 2, -0.95,
        r"(b) returning to the past requires $v > c$",
        ha="center", va="top", fontsize=9)

plt.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "worldline_figure.pdf")
plt.savefig(out, bbox_inches="tight")
print(f"Wrote {out}")
