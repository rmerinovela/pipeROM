"""
Equivalent SDOF vs 3D NLTHA: peak displacements and displaced shapes.

Used by the <model>_biron.ipynb notebooks of this folder. Data read:
  3D NLTHA     Results/<model>_biron_Disp<X|Y>.txt        peak displacement over the braced nodes
                                                          (rows: Names.txt records, cols: IM1-IM10;
                                                          0 = analysis not available)
               Results/<model>_bironDispShape<X|Y>.npy    normalized displaced shapes (record, node, IM)
  SDOF NLTHA   Pushover_SDOF/Results<tag>_SDOF/<tag>_NLTHA_peak_displacements.csv
  2D pushover  Pushover2D/pushover_results_<tag>.txt, pushover_layout_<tag>.json, and nltha_lines /
               x_ref of Pushover2D/CS_Lumped_Iter_<tag>_cont_PO.py (NLTHA points on the analysed line)
  3D geometry  inputs/models3d/<model>.json               node coordinates, elements, DispShape nodes
"""
import os
import io
import ast
import csv
import sys
import json
import contextlib
from collections import deque

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle

RESULTS_DIR    = os.path.dirname(os.path.abspath(__file__))
PAPER_DIR      = os.path.dirname(RESULTS_DIR)
SDOF_DIR       = os.path.join(PAPER_DIR, "Pushover_SDOF")
PUSHOVER2D_DIR = os.path.join(PAPER_DIR, "Pushover2D")
MODELS3D_DIR   = os.path.join(PAPER_DIR, "..", "inputs", "models3d")

N_RECORDS = 44
IMS       = np.arange(1, 11)          # intensity levels of the 3D NLTHA

# Collapse (3D_models/*_biron.py): a braced node beyond COLLAPSE_DISP stops the run, whose two records
# are stored with COLLAPSE_PEAK_3D. The same rule is applied to the SDOF (peak_comparison), with
# collapsed records stored as SDOF_COLLAPSE_PEAK.
COLLAPSE_DISP      = 60.0     # mm
COLLAPSE_PEAK_3D   = 100.0    # mm
SDOF_COLLAPSE_PEAK = 100.0    # mm

sys.path.insert(0, SDOF_DIR)
from sdof_from_2d import load_sdof_params   # noqa: E402

C_3D, C_SDOF, C_STATIC = "k", "#d1495b", "C0"


# ----------------------------------------------------------------------------------------------
# Peak displacements
# ----------------------------------------------------------------------------------------------
def load_3d_peaks(model, direction):
    """3D NLTHA peak displacements (mm), shape (44 records, 10 IMs); NaN where not available."""
    d = np.loadtxt(os.path.join(RESULTS_DIR, f"{model}_biron_Disp{direction.upper()}.txt"))
    return np.where(d > 0, d, np.nan)


def load_sdof_peaks(model, direction):
    """SDOF NLTHA peak displacements (mm), shape (44, 10); NaN for IMs not yet run or records
    that did not reach the end of the motion without collapsing. Collapsed records (stopped when a
    braced node exceeded the collapse displacement) keep their peak, beyond it by construction.
    Rows follow the records' order (Names.txt)."""
    tag  = model + direction.lower()
    path = os.path.join(SDOF_DIR, f"Results{tag}_SDOF", f"{tag}_NLTHA_peak_displacements.csv")
    out  = np.full((N_RECORDS, len(IMS)), np.nan)
    if not os.path.exists(path):
        return out
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    count = {}
    for r in rows:
        im = int(r["IM"])
        k  = count.get(im, 0)
        count[im] = k + 1
        if im <= len(IMS) and (int(r["completed"]) or int(r.get("collapsed") or 0)):
            out[k, im - 1] = float(r["peak_disp_mm"])
    return out


def sdof_params(model, direction, dc_target=12.0):
    """Equivalent SDOF parameters (load_sdof_params) without its printout."""
    with contextlib.redirect_stdout(io.StringIO()):
        return load_sdof_params(model + direction.lower(), dc_target=dc_target)


def sdof_to_braced_factor(model, direction, dc_target=12.0):
    """Gamma * max(phi) over the spring DOFs (stiff hangers + orthogonals): converts the SDOF
    displacement to the largest displacement at the braced nodes, the quantity of Disp<X|Y>.txt."""
    p = sdof_params(model, direction, dc_target)
    return p["Gamma"] * max(phi for phi, _ in p["trans"] + p["long"])


def peak_comparison(model, factor="gamma_phi", dc_target=12.0, sdof_collapse=True):
    """
    Peaks of both directions: {dir: dict(d3=..., sdof=..., factor=..., sdof_collapsed=...)}; sdof
    already multiplied by the factor. factor = "gamma_phi" (sdof_to_braced_factor), a number, or a
    dict {"X": .., "Y": ..}.

    sdof_collapse: apply the 3D collapse rule to the SDOF. A 3D run applies record r in x and record
    r^1 (its pair, 2g <-> 2g+1) in y, and stops both when a braced node exceeds COLLAPSE_DISP in either
    direction. Likewise, an SDOF pair (X row r, Y row r^1) collapses when either braced-node peak
    exceeds COLLAPSE_DISP, and both are then stored as SDOF_COLLAPSE_PEAK.
    """
    out = {}
    for d in "XY":
        f = factor[d] if isinstance(factor, dict) else factor
        if f == "gamma_phi":
            f = sdof_to_braced_factor(model, d, dc_target)
        out[d] = dict(d3=load_3d_peaks(model, d), sdof=float(f) * load_sdof_peaks(model, d),
                      factor=float(f))
    pair = np.arange(N_RECORDS) ^ 1
    over = {d: out[d]["sdof"] > COLLAPSE_DISP for d in "XY"}    # NaN compares False
    coll = {"X": over["X"] | over["Y"][pair], "Y": over["Y"] | over["X"][pair]}
    for d in "XY":
        out[d]["sdof_collapsed"] = coll[d] if sdof_collapse else np.zeros_like(coll[d])
        if sdof_collapse:
            out[d]["sdof"] = np.where(coll[d], SDOF_COLLAPSE_PEAK, out[d]["sdof"])
    return out


def collapsed_3d_median(comp, d):
    """True for the IMs where at least half of the available 3D records collapsed (COLLAPSE_PEAK_3D),
    i.e. the 3D median is a collapse or, with exactly half, the mean of a collapse and a real peak."""
    a = comp[d]["d3"]
    n = np.sum(~np.isnan(a), axis=0)
    return (n > 0) & (np.sum(a >= COLLAPSE_PEAK_3D - 1e-6, axis=0) >= 0.5 * n)


def _median(a):
    """Median over records of each IM (NaN where an IM has no data)."""
    m = np.full(a.shape[1], np.nan)
    for j in range(a.shape[1]):
        v = a[~np.isnan(a[:, j]), j]
        if len(v):
            m[j] = np.median(v)
    return m


def _pct(a, q):
    p = np.full(a.shape[1], np.nan)
    for j in range(a.shape[1]):
        v = a[~np.isnan(a[:, j]), j]
        if len(v):
            p[j] = np.percentile(v, q)
    return p


def peak_table(comp):
    """Medians per IM and direction, with the SDOF/3D ratio (pandas DataFrame)."""
    import pandas as pd
    cols = {}
    for d, c in comp.items():
        m3, ms = _median(c["d3"]), _median(c["sdof"])
        cols[(d, "3D (mm)")]     = m3
        cols[(d, "SDOF (mm)")]   = ms
        cols[(d, "SDOF / 3D")]   = ms / m3
        cols[(d, "n 3D / SDOF")] = [f"{np.sum(~np.isnan(c['d3'][:, j]))} / {np.sum(~np.isnan(c['sdof'][:, j]))}"
                                    for j in range(len(IMS))]
    df = pd.DataFrame(cols, index=pd.Index(IMS, name="IM"))
    return df.round(3)


def plot_peak_medians(model, comp):
    """Median (markers) and 16th-84th percentiles (bars) of the peak displacement vs IM."""
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (d, c) in zip(axs, comp.items()):
        for a, color, mk, dx, label in [(c["d3"], C_3D, "o", -0.08, "3D NLTHA"),
                                        (c["sdof"], C_SDOF, "s", 0.08,
                                         rf"SDOF NLTHA $\times$ {c['factor']:.3f}")]:
            m, lo, hi = _median(a), _pct(a, 16), _pct(a, 84)
            ax.errorbar(IMS + dx, m, yerr=[m - lo, hi - m], fmt=mk, color=color, ms=6,
                        capsize=3, lw=1, label=label)
        ax.set_title(f"{model} – direction {d}")
        ax.set_xlabel("Intensity level", fontsize=13)
        ax.set_ylabel("Peak displacement at braced nodes (mm)", fontsize=12)
        ax.set_xticks(IMS)
        # Medians in view; the upper percentiles of collapsing records may extend beyond
        top = np.nanmax(np.concatenate([_median(c["d3"]), _median(c["sdof"]), _pct(c["sdof"], 84)]))
        ax.set_ylim(0, 1.4 * top)
        ax.legend(fontsize=10, loc="upper left")
    fig.tight_layout()
    return fig


def plot_peak_ratio(model, comp):
    """Ratio of the medians, SDOF / 3D, vs IM with the ±10/20/30% bands."""
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for b, ls, color in [(0.1, "--", "0.6"), (0.2, "--", "k"), (0.3, "-.", "k")]:
        ax.axhline(1 + b, color=color, ls=ls, lw=1)
        ax.axhline(1 - b, color=color, ls=ls, lw=1)
    ax.axhline(1, color="k", lw=1)
    for d, mk in zip(comp, "os"):
        c = comp[d]
        ax.plot(IMS, _median(c["sdof"]) / _median(c["d3"]), mk, ms=7, label=f"Direction {d}")
    ax.set_xlim(0, IMS[-1] + 0.5)
    ax.set_ylim(0, 2)
    ax.set_xticks(IMS)
    ax.set_xlabel("Intensity level", fontsize=14)
    ax.set_ylabel(r"$\tilde{D}_{SDOF}\,/\,\tilde{D}_{3D}$", fontsize=16)
    ax.set_title(f"{model}: ratio of median peak displacements")
    ax.legend(fontsize=11)
    fig.tight_layout()
    return fig


def plot_peak_records(model, comp):
    """Record by record: SDOF vs 3D peak displacement, coloured by IM."""
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
    cmap = plt.get_cmap("viridis", len(IMS))
    for ax, (d, c) in zip(axs, comp.items()):
        ok = ~np.isnan(c["d3"]) & ~np.isnan(c["sdof"])
        if not ok.any():
            ax.text(0.5, 0.5, "no SDOF results yet", ha="center", transform=ax.transAxes)
            continue
        ims = np.broadcast_to(IMS, c["d3"].shape)
        sc = ax.scatter(c["d3"][ok], c["sdof"][ok], c=ims[ok], cmap=cmap, vmin=0.5,
                        vmax=IMS[-1] + 0.5, s=14, alpha=0.8)
        lo = 0.8 * min(c["d3"][ok].min(), c["sdof"][ok].min())
        hi = 1.25 * max(c["d3"][ok].max(), c["sdof"][ok].max())
        xs = np.array([lo, hi])
        ax.plot(xs, xs, "k", lw=1)
        ax.plot(xs, 1.2 * xs, "k--", lw=0.8)
        ax.plot(xs, xs / 1.2, "k--", lw=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xlabel("3D NLTHA peak (mm)", fontsize=12)
        ax.set_ylabel(rf"SDOF NLTHA peak $\times$ {c['factor']:.3f} (mm)", fontsize=12)
        ax.set_title(f"{model} – direction {d} ({ok.sum()} record–IM pairs)")
        fig.colorbar(sc, ax=ax, label="Intensity level", ticks=IMS)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------------------------
# Displaced shapes
# ----------------------------------------------------------------------------------------------
def load_geometry(model):
    """
    Plan geometry of the DispShape nodes: dict with
      xy    : (n, 2) coordinates (m), in the order of the DispShape arrays
      edges : pairs of DispShape node indices joined by pipe (through any non-recorded nodes)
      braced: {"X": [...], "Y": [...]} DispShape indices of the braced nodes of each direction
    """
    with open(os.path.join(MODELS3D_DIR, f"{model}.json")) as f:
        d = json.load(f)
    coords = {int(a[0]): np.array(a[1:4], float) for c, a in d["commands"] if c == "node"}
    nodes  = [int(n) for n in d["node_groups"]["X"]]
    z_pipe = coords[nodes[0]][2]

    # Pipe network: elements between nodes at the pipe level (pipes and pipe joints)
    adj = {}
    for c, a in d["commands"]:
        if c != "element":
            continue
        i, j = int(a[2]), int(a[3])
        if i in coords and j in coords and coords[i][2] == z_pipe and coords[j][2] == z_pipe:
            adj.setdefault(i, set()).add(j)
            adj.setdefault(j, set()).add(i)

    # Contract to the recorded nodes: walk from each one until another recorded node is reached
    index = {n: k for k, n in enumerate(nodes)}
    edges = set()
    for n in nodes:
        seen, todo = {n}, deque(adj.get(n, ()))
        while todo:
            q = todo.popleft()
            if q in seen:
                continue
            seen.add(q)
            if q in index:
                edges.add(tuple(sorted((index[n], index[q]))))
            else:
                todo.extend(adj.get(q, ()))
    xy = np.array([coords[n][:2] for n in nodes]) / 1000
    # Braced nodes (transverse restraints) of each direction, as DispShape indices
    braced = {dr: [index[n] for n in d["node_groups"][dr + "t"] if n in index] for dr in "XY"}
    return dict(xy=xy, edges=sorted(edges), nodes=nodes, braced=braced)


def _driver_inputs(tag):
    """nltha_lines and x_ref (mm, None if not given) of the 2D pushover driver script."""
    path = os.path.join(PUSHOVER2D_DIR, f"CS_Lumped_Iter_{tag}_cont_PO.py")
    with open(path) as f:
        tree = ast.parse(f.read())
    lines, x_ref = None, None
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "nltha_lines" for t in node.targets):
            lines = ast.literal_eval(node.value)
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "plot_shape_vs_nltha":
            for kw in node.keywords:
                if kw.arg == "x_ref":
                    x_ref = ast.literal_eval(kw.value)
    if lines is None:
        raise ValueError(f"No nltha_lines in {path}")
    return lines, x_ref


def static_shape(model, direction, geom, dc_target=12.0):
    """
    Displaced shape of the equivalent static procedure (2D pushover at the step closest to
    dc_target) at the DispShape nodes.

    The 2D model is a line (the NLTHA points of nltha_lines, with their 2D coordinate). Every
    other node takes the 2D coordinate of the node it hangs from: the same one along pipes
    parallel to the excitation (rigid axially), shifted by the distance along pipes parallel
    to the 2D line (branches not modelled deform as the analysed one).

    Returns dict(phi=values per node, lines=[(label, idx)], i_ref=normalization node, dc=step).
    """
    tag = model + direction.lower()
    res = np.loadtxt(os.path.join(PUSHOVER2D_DIR, f"pushover_results_{tag}.txt"))
    with open(os.path.join(PUSHOVER2D_DIR, f"pushover_layout_{tag}.json")) as f:
        layout = json.load(f)
    x_d = np.asarray(layout["x_d"], float)
    n   = len(x_d)
    row = res[np.argmin(np.abs(res[:, 0] - dc_target))]
    d_norm = np.empty(n)
    d_norm[layout["column_order"]] = row[6:6 + n]

    # 2D shape as a function of the 2D coordinate (DOFs at the same x averaged)
    xs = np.unique(x_d)
    ds = np.array([d_norm[x_d == x].mean() for x in xs])
    shape2d = lambda x: np.interp(x, xs, ds)

    lines, x_ref = _driver_inputs(tag)
    xy  = geom["xy"] * 1000
    x2d = np.full(len(xy), np.nan)
    for _, idx, xl in lines:
        x2d[idx] = xl

    # Direction of the 2D line in plan (mm of plan per mm of 2D coordinate)
    _, idx, xl = lines[0]
    a, b = np.argmin(xl), np.argmax(xl)
    u = (xy[idx[b]] - xy[idx[a]]) / (xl[b] - xl[a])
    e = np.array([1.0, 0.0]) if direction.upper() == "X" else np.array([0.0, 1.0])

    nb = {}
    for i, j in geom["edges"]:
        nb.setdefault(i, []).append(j)
        nb.setdefault(j, []).append(i)
    todo = deque(np.flatnonzero(~np.isnan(x2d)))
    while todo:
        p = todo.popleft()
        for q in nb.get(p, []):
            if np.isnan(x2d[q]):
                dxy = xy[q] - xy[p]
                parallel = abs(dxy @ e) > 0.9 * np.linalg.norm(dxy)
                x2d[q] = x2d[p] if parallel else x2d[p] + dxy @ u / (u @ u)
                todo.append(q)

    phi = shape2d(x2d)

    # Normalization node: the analysed-line node closest to x_ref (default: last 2D DOF)
    if x_ref is None:
        x_ref = x_d[-1]
    on_line = np.concatenate([idx for _, idx, _ in lines])
    i_ref   = on_line[np.argmin(np.abs(x2d[on_line] - x_ref))]
    phi     = phi / phi[i_ref]
    return dict(phi=phi, lines=[(lab, np.asarray(idx)) for lab, idx, _ in lines], i_ref=i_ref,
                dc=row[0], x2d=x2d)


def nltha_shape(model, direction, im, i_ref):
    """Mean and standard deviation of the 3D NLTHA displaced shapes at intensity level im
    (records with all-zero shapes excluded), normalized so that the mean is 1 at node i_ref."""
    s = np.load(os.path.join(RESULTS_DIR, f"{model}_bironDispShape{direction.upper()}.npy"))[:, :, im - 1]
    s = s[np.any(s > 0, axis=1)]
    mu, sd = s.mean(axis=0), s.std(axis=0)
    return mu / mu[i_ref], sd / mu[i_ref], len(s)


def _draw(ax, geom, e, amp, mu, sd, phi):
    xy, edges = geom["xy"], geom["edges"]
    segs = lambda v: [np.array([xy[i] + amp * v[i] * e, xy[j] + amp * v[j] * e]) for i, j in edges]

    band = [np.array([xy[i] + amp * (mu[i] - 2 * sd[i]) * e, xy[j] + amp * (mu[j] - 2 * sd[j]) * e,
                      xy[j] + amp * (mu[j] + 2 * sd[j]) * e, xy[i] + amp * (mu[i] + 2 * sd[i]) * e])
            for i, j in edges]
    ax.add_collection(PolyCollection(band, facecolor="0.5", edgecolor="0.5", lw=0.3, alpha=0.45,
                                     label=r"NLTHA $\mu \pm 2\sigma$"))
    for k, (p, q) in enumerate(segs(np.zeros(len(xy)))):
        ax.plot(*np.array([p, q]).T, color="k", lw=2, label="Original position" if k == 0 else None)
    for k, (p, q) in enumerate(segs(mu)):
        ax.plot(*np.array([p, q]).T, color="0.35", ls="-.", lw=1.6,
                label=r"Displaced shape NLTHA ($\mu$)" if k == 0 else None)
    for k, (p, q) in enumerate(segs(phi)):
        ax.plot(*np.array([p, q]).T, color=C_STATIC, ls="--", lw=1.6,
                label="Displaced shape static" if k == 0 else None)


def plot_displaced_shape(model, direction, im, amp=None, amp_zoom=None, dc_target=12.0):
    """
    Plan view of the NLTHA mean displaced shape (±2σ band) and of the equivalent static shape,
    both normalized to 1 at the same node (marked). Top: whole system, displacements amplified
    by amp (m per unit normalized displacement). Bottom: each analysed line of the 2D model,
    zoomed (amplification amp_zoom, default amp).
    """
    direction = direction.upper()
    geom = load_geometry(model)
    st   = static_shape(model, direction, geom, dc_target)
    mu, sd, n_rec = nltha_shape(model, direction, im, st["i_ref"])
    phi  = st["phi"]
    xy   = geom["xy"]
    e    = np.array([1.0, 0.0]) if direction == "X" else np.array([0.0, 1.0])

    extent = np.ptp(xy, axis=0).max()
    if amp is None:
        amp = 0.08 * extent / max(np.nanmax(np.abs(mu)), np.nanmax(np.abs(phi)))
    amp_zoom = amp if amp_zoom is None else amp_zoom

    n_lines = len(st["lines"])
    fig = plt.figure(figsize=(11, 10.5))
    gs  = fig.add_gridspec(2, n_lines, height_ratios=[1.35, 1])
    ax  = fig.add_subplot(gs[0, :])

    _draw(ax, geom, e, amp, mu, sd, phi)
    ax.plot(*xy[st["i_ref"]], "o", mfc="w", mec="k", ms=7, zorder=5, label="Normalization node")
    ax.set_aspect("equal")
    ax.autoscale_view()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    pad = 0.06 * extent
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_xlabel(r"$x$  (m)", fontsize=13)
    ax.set_ylabel(r"$y$  (m)", fontsize=13)
    ax.set_title(f"{model} – direction {direction}: NLTHA at IM{im} ({n_rec} records) vs static "
                 rf"($\Delta_c$ = {st['dc']:.1f} mm); displacements $\times${amp:.3g} m")
    ax.legend(fontsize=10, loc="upper left", bbox_to_anchor=(1.01, 1))

    for k, (label, idx) in enumerate(st["lines"]):
        az = fig.add_subplot(gs[1, k])
        _draw(az, geom, e, amp_zoom, mu, sd, phi)
        pts = np.vstack([xy[idx] + amp_zoom * np.outer(v, e)
                         for v in (np.zeros(len(idx)), mu[idx] - 2 * sd[idx], mu[idx] + 2 * sd[idx], phi[idx])])
        lo, hi = pts.min(axis=0), pts.max(axis=0)
        margin = np.maximum(0.08 * (hi - lo), 0.02 * extent)
        lo, hi = lo - margin, hi + margin
        az.set_xlim(lo[0], hi[0])
        az.set_ylim(lo[1], hi[1])
        az.set_xlabel(r"$x$  (m)", fontsize=12)
        az.set_ylabel(r"$y$  (m)", fontsize=12)
        tag = f"({chr(ord('a') + k)}) " if n_lines > 1 else ""
        az.set_title(f"{tag}Zoom: analysed line {label}".strip() +
                     (rf", displacements $\times${amp_zoom:.3g} m" if amp_zoom != amp else ""), fontsize=11)
        if amp_zoom == amp:
            ax.add_patch(Rectangle(lo, *(hi - lo), fill=False, ec="C3", lw=1.2, ls=":"))
            ax.annotate(tag.strip() or "zoom", hi, color="C3", fontsize=10)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------------------------
# Pushover curves
# ----------------------------------------------------------------------------------------------
C_2D = "C0"


def load_pushovers(model, direction, dc_target=12.0):
    """
    Pushover curves of one direction as (u_ref, Vb) with u_ref = displacement (mm) at the reference
    node of the 2D model (where its shape is normalized) and Vb = base shear (N):
      2d   Pushover2D/pushover_results_<tag>.txt: u_ref = Gamma * u_SDOF (= Delta_c / max d_norm)
      sdof Pushover_SDOF/Results<tag>_SDOF/DispC.out, VbaseC.out: u_ref = Gamma * u (phi_ref = 1),
           Vb = - support reaction
      3d   Results/<model>_biron_pushover<X|Y>.txt (3D_models/pushover_3d.py): displacement of the
           3D node at the reference DOF, under a mass-proportional load; None if not run yet
    and the 2D step the SDOF was built from (calib = (u_ref, Vb, Delta_c)).
    """
    tag = model + direction.lower()
    res = np.loadtxt(os.path.join(PUSHOVER2D_DIR, f"pushover_results_{tag}.txt"))
    out = {"2d": (np.r_[0.0, res[:, 1] * res[:, 5]], np.r_[0.0, res[:, 3]])}
    i = np.argmin(np.abs(res[:, 0] - dc_target))
    out["calib"] = (res[i, 1] * res[i, 5], res[i, 3], res[i, 0])

    d = os.path.join(SDOF_DIR, f"Results{tag}_SDOF")
    gamma = sdof_params(model, direction, dc_target)["Gamma"]
    u, r = np.loadtxt(os.path.join(d, "DispC.out")), np.loadtxt(os.path.join(d, "VbaseC.out"))
    n = min(len(u), len(r))
    out["sdof"] = (np.r_[0.0, gamma * u[:n]], np.r_[0.0, -r[:n]])

    p = os.path.join(RESULTS_DIR, f"{model}_biron_pushover{direction.upper()}.txt")
    out["3d"] = tuple(np.loadtxt(p).T) if os.path.exists(p) else None
    return out


def plot_pushover_comparison(model, dc_target=12.0, u_max=None):
    """Base shear vs displacement at the reference node of the 2D model: 3D pushover (mass-proportional
    load), 2D adaptive pushover and SDOF pushover, both directions."""
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
    for ax, d in zip(axs, "XY"):
        c = load_pushovers(model, d, dc_target)
        if c["3d"] is not None:
            ax.plot(c["3d"][0], c["3d"][1] / 1e3, color=C_3D, lw=1.8, label="3D pushover")
        else:
            ax.text(0.5, 0.5, "3D pushover not run yet\n(3D_models/pushover_3d.py)", ha="center",
                    transform=ax.transAxes, color="0.4")
        ax.plot(c["2d"][0], c["2d"][1] / 1e3, "o--", color=C_2D, ms=3.5, lw=1.4, label="2D adaptive pushover")
        ax.plot(c["sdof"][0], c["sdof"][1] / 1e3, color=C_SDOF, lw=1.4, label=r"SDOF pushover ($\Gamma u$)")
        u, v, dc = c["calib"]
        ax.plot(u, v / 1e3, "*", color=C_SDOF, mec="k", ms=13, zorder=5,
                label=rf"SDOF calibration ($\Delta_c$ = {dc:.1f} mm)")
        ax.set_xlim(0, u_max or max(c["2d"][0].max(), 0 if c["3d"] is None else c["3d"][0].max()))
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Displacement at the reference node (mm)", fontsize=12)
        ax.set_ylabel("Base shear (kN)", fontsize=12)
        ax.set_title(f"{model} – direction {d}")
        ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    return fig
