"""Plan view of archetype 2 (y direction) divided into tributary areas (worked example, Figure 1).

Run from the repository root: python examples/plot_tributary_areas.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

OUT = Path(__file__).parent / "figures" / "archetype2_y_tributary_areas.png"

L = 36.0                                   # main (transverse) line, m
HANGERS = [1.0 + 3.0 * i for i in range(12)]
RESTRAINTS = [1.0, 13.0, 22.0, 34.0]       # transverse restraints (mask of M02y.yaml)
M_PER_M = 0.1248                           # t/m, steel + water x 1.35

C_HANGER = ("#c6dbef", "#9ecae1")
C_RESTR = "#fdae6b"
C_LONG = "#a1d99b"
C_LINE = "#404040"

# Tributary limits: half-way between neighbouring nodes; the junction keeps 35-36 m
nodes = HANGERS + [L]
edges = [0.5] + [(a + b) / 2 for a, b in zip(nodes[:-1], nodes[1:])] + [L]

fig, ax = plt.subplots(figsize=(12, 4.2))
H = 0.9                                    # half-height of the tributary bands (drawing units)

# Free end, outside the procedure
ax.add_patch(Rectangle((0, -H), 0.5, 2 * H, facecolor="none", edgecolor="0.6", hatch="////", lw=0.6))

for i, (x, a, b) in enumerate(zip(nodes, edges[:-1], edges[1:])):
    if x == L:
        color = C_LONG
    elif x in RESTRAINTS:
        color = C_RESTR
    else:
        color = C_HANGER[i % 2]
    ax.add_patch(Rectangle((a, -H), b - a, 2 * H, facecolor=color, edgecolor="white", lw=1.0))
    if x == L:
        continue
    ax.text((a + b) / 2, -H - 0.35, f"{(b - a) * M_PER_M:.3f}", ha="center", va="top", fontsize=8)
    ax.text((a + b) / 2, -H - 1.05, f"{b - a:g} m", ha="center", va="top", fontsize=7, color="0.35")

ax.text(-0.4, -H - 0.35, "mass (t)", ha="right", va="top", fontsize=8)
ax.text(-0.4, -H - 1.05, "tributary length", ha="right", va="top", fontsize=7, color="0.35")

# Main line, hangers and transverse restraints
ax.plot([0, L], [0, 0], color=C_LINE, lw=3, solid_capstyle="butt", zorder=3)
for x in HANGERS:
    if x in RESTRAINTS:
        ax.add_patch(Polygon([[x, 0.1], [x - 0.45, 0.75], [x + 0.45, 0.75]], closed=True,
                             facecolor="#d94801", edgecolor="black", lw=0.6, zorder=4))
        ax.add_patch(Polygon([[x, -0.1], [x - 0.45, -0.75], [x + 0.45, -0.75]], closed=True,
                             facecolor="#d94801", edgecolor="black", lw=0.6, zorder=4))
    ax.plot(x, 0, "o", ms=5, mfc="white", mec="black", mew=0.8, zorder=5)
    ax.text(x, H + 0.15, f"{x:g}", ha="center", va="bottom", fontsize=7, color="0.35")
ax.text(L + 0.7, H + 0.15, "36", ha="left", va="bottom", fontsize=7, color="0.35")
ax.text(-0.4, H + 0.15, "x (m)", ha="right", va="bottom", fontsize=7, color="0.35")

# Longitudinal line (branch), drawn shortened
BR = 7.0
ax.add_patch(Rectangle((L - 0.5, H), 1.0, BR - H, facecolor=C_LONG, edgecolor="none", alpha=0.6))
ax.plot([L, L], [0, BR], color=C_LINE, lw=3, solid_capstyle="butt", zorder=3)
for yb in (3.0, 5.6):                      # positions are schematic: only the number of restraints is used
    ax.plot([L - 0.55, L + 0.55], [yb, yb], color="#238b45", lw=3, zorder=4)
for yb in (4.15, 4.45):                    # break marks: the branch is 18 m long
    ax.plot([L - 0.5, L + 0.5], [yb - 0.15, yb + 0.15], color="black", lw=0.8, zorder=6)
ax.plot(L, 0, "s", ms=7, mfc="#238b45", mec="black", mew=0.8, zorder=5)
ax.annotate("longitudinal line (18 m, 2 longitudinal restraints)\n"
            r"$m_{L,1}$ = 2.246 t on spring $k_{L,1}$",
            xy=(L + 0.6, 5.6), xytext=(L + 1.6, 6.2), fontsize=8, va="center",
            arrowprops=dict(arrowstyle="-", lw=0.6))
ax.annotate("junction (1 m):\n" r"$m_{T|L,1}$ = 0.125 t", xy=(35.5, -H), xytext=(L + 1.6, -2.0),
            fontsize=8, va="center", arrowprops=dict(arrowstyle="-", lw=0.6))
ax.annotate("free end (0.06 t),\nnot part of the procedure", xy=(0.25, -H), xytext=(-6.3, -3.9),
            fontsize=7, va="center", color="0.35", arrowprops=dict(arrowstyle="-", lw=0.6, color="0.5"))

# Loading direction
ax.annotate("", xy=(1.0, 5.8), xytext=(1.0, 3.4), arrowprops=dict(arrowstyle="-|>", lw=1.4))
ax.text(1.4, 4.6, "y (loading)", fontsize=8, va="center")
ax.annotate("", xy=(3.4, 3.4), xytext=(1.0, 3.4), arrowprops=dict(arrowstyle="-|>", lw=1.0, color="0.4"))
ax.text(3.7, 3.4, "x", fontsize=8, va="center", color="0.4")

# Legend
items = [(C_HANGER[1], "hanger"), (C_RESTR, "hanger with transverse restraint"),
         (C_LONG, "longitudinal line and junction")]
for k, (c, label) in enumerate(items):
    x0 = (9.0, 16.0, 28.5)[k]
    ax.add_patch(Rectangle((x0, -4.4), 0.9, 0.5, facecolor=c, edgecolor="none"))
    ax.text(x0 + 1.2, -4.15, label, fontsize=8, va="center")

ax.set_xlim(-6.5, L + 13)
ax.set_ylim(-4.8, 7.4)
ax.axis("off")
OUT.parent.mkdir(exist_ok=True)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print(OUT)
