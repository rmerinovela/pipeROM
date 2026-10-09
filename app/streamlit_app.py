"""Interactive front end for the reduced-order model of suspended piping systems.

The user describes the whole system in plan (piperom.layout): a main line along x and identical branches on
each side, with restraints in both directions. The app derives the equivalent system of each loading
direction and runs every analysis for both.

    streamlit run app/streamlit_app.py
"""

from __future__ import annotations

import hashlib
import io
import multiprocessing as mp
import sys
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from piperom.inputs import AnalysisSettings, InputError, PipingSystem, dump_yaml  # noqa: E402
from piperom.jobs import pushover_job, sdof_job, sdof_time_history_job, verification_job  # noqa: E402
from piperom.layout import Layout, layout_from_dict  # noqa: E402
from piperom.sdof import sdof_delta_c  # noqa: E402
from piperom.motions import load_motion_sets, select_runs  # noqa: E402
from piperom.timehistory import SDOFModel, support_demands  # noqa: E402
from piperom.verification3d import braced_peaks, list_models3d, load_model3d  # noqa: E402
from piperom.pushover import CURVE_COLUMNS, SHAPE_COLUMNS, rows_to_csv  # noqa: E402
from piperom.trapeze import (PARAMETER_NAMES, Pinching4, load_trapeze, parse_trapeze_csv,  # noqa: E402
                             trilinear_from_pinching4)

INPUT_DIRS = {"Archetype": REPO / "inputs" / "layouts", "Example": REPO / "automatized_workflow" / "examples"}
KINDS = {"transverse": "Transverse (C-TPS-T)", "longitudinal": "Longitudinal (C-TPS-L)"}
PREFIX = {"transverse": "T", "longitudinal": "L"}
DIRS = ("x", "y")
SIDE_KEY = {"+y": "p", "-y": "n"}
PIPE_KEYS = ("outer_diameter", "inner_diameter", "elastic_modulus", "shear_modulus", "density", "fluid_density",
             "mass_factor")
C_MAIN, C_BRACE, C_BRACE_X, C_BRANCH, C_MUTED = "#1f77b4", "#d62728", "#9467bd", "#2ca02c", "#9e9e9e"

st.set_page_config(page_title="Suspended piping ROM", layout="wide")


# ============================================================================ worker process
@st.cache_resource
def executor() -> ProcessPoolExecutor:
    # OpenSees state is global per process: analyses run in worker processes, one at a time each.
    return ProcessPoolExecutor(max_workers=2, mp_context=mp.get_context("spawn"))


@st.cache_resource
def manager():
    return mp.get_context("spawn").Manager()


# ============================================================================ state
def _hanger_state(prefix: str, spec: dict, positions: np.ndarray) -> None:
    ss = st.session_state
    if spec.get("positions") is None and spec.get("first") is not None:
        ss[f"{prefix}_hmode"] = "Regular grid"
        ss[f"{prefix}_hfirst"] = float(spec["first"])
        ss[f"{prefix}_hspacing"] = float(spec["spacing"])
    else:
        ss[f"{prefix}_hmode"] = "Explicit positions"
        ss[f"{prefix}_hfirst"], ss[f"{prefix}_hspacing"] = 1000.0, 3000.0
    ss[f"{prefix}_hclear"] = float(spec.get("end_clearance", 1000.0))
    ss[f"{prefix}_hpos"] = ", ".join(f"{x:g}" for x in positions)


def load_into_state(layout: Layout, trapezes: dict[str, Pinching4] | None = None) -> None:
    ss = st.session_state
    raw = layout.raw
    ss.v = ss.get("v", 0) + 1
    ss.name, ss.description = layout.name, layout.description
    pipe = PipingSystem.from_dict({"main_line": {"length": 1000, "n_pipes": 1}, "braces": {"count": 0},
                                   "pipe": layout.pipe}).pipe
    for k in PIPE_KEYS:
        ss[f"pipe_{k}"] = getattr(pipe, k)
    ss.branch_participation = 1.0
    m = layout.main
    ss.m_length, ss.m_npipes = float(m.length), int(m.n_pipes)
    _hanger_state("m", raw["main_line"].get("hangers") or {}, m.hangers)
    ss.m_brace_t = {float(x) for x in m.hangers[m.transverse == 1]}
    ss.m_brace_l = {float(x) for x in m.hangers[m.longitudinal == 1]}
    first = layout.branches[0].line
    ss.b_length, ss.b_npipes = float(first.length), int(first.n_pipes)
    raw_br = raw.get("branches") or []
    for side, s in SIDE_KEY.items():
        bs = [b for b in layout.branches if b.side == side]
        ss[f"{s}_tees"] = ", ".join(f"{b.x:g}" for b in bs)
        line = bs[0].line if bs else first
        spec = next((rb.get("hangers") for rb in raw_br if rb.get("side") == side and rb.get("hangers")), None)
        _hanger_state(s, spec or raw["main_line"].get("hangers") or {}, line.hangers)
        ss[f"{s}_brace_t"] = {float(x) for x in line.hangers[line.transverse == 1]} if bs else set()
        ss[f"{s}_brace_l"] = {float(x) for x in line.hangers[line.longitudinal == 1]} if bs else set()
    trapezes = trapezes or {}
    ss.trapeze_base = {k: trapezes.get(k) or load_trapeze(layout.trapezes.get(k, "default"), k) for k in KINDS}
    ss.trapeze_source = {k: "default file" if layout.trapezes.get(k, "default") == "default"
                         else Path(str(layout.trapezes[k])).name for k in KINDS}


def load_settings_into_state(s: AnalysisSettings) -> None:
    ss = st.session_state
    ss.s_explicit = s.delta_c_values is not None
    ss.s_values = ", ".join(f"{v:g}" for v in (s.delta_c_values or []))
    ss.s_start, ss.s_stop, ss.s_n = s.delta_c_start, s.delta_c_stop, s.n_steps
    ss.s_warm = s.warm_start
    ss.s_maxit, ss.s_tol = s.max_iterations, s.tolerance
    ss.s_test, ss.s_stol, ss.s_smaxit = s.solver_test, s.solver_tolerance, s.solver_max_iterations
    ss.s_sdof_dc, ss.s_sdof_grid = s.sdof_delta_c, s.sdof_on_pushover_grid
    ss.s_split = s.branch_split
    ss.s_raw = {k: s.to_dict()[k] for k in ("sdof", "motions", "sdof_time_history", "verification_3d")}
    th, v3, mo = s.sdof_time_history, s.verification_3d, s.motions
    ss.th_xi, ss.th_test, ss.th_tol, ss.th_maxit = th["damping_ratio"], th["test"], th["tolerance"], th["max_iterations"]
    ss.v3_xi, ss.v3_test, ss.v3_tol, ss.v3_maxit = v3["damping_ratio"], v3["test"], v3["tolerance"], v3["max_iterations"]
    ss.v3_fb_test, ss.v3_fb_tol = v3["fallback_test"], v3["fallback_tolerance"]
    ss.m_set, ss.m_levels, ss.m_floor = mo.get("set"), list(mo.get("levels") or []), int(mo.get("floor", 4))


def load_layout_file(path: Path) -> Layout:
    return layout_from_dict(yaml.safe_load(path.read_text()) or {}, base_dir=path.parent)


if "v" not in st.session_state:
    load_into_state(load_layout_file(INPUT_DIRS["Archetype"] / "M02.yaml"))
    load_settings_into_state(AnalysisSettings.from_dict())

# Streamlit drops the state of widgets that are not rendered in a run (e.g. grid inputs while explicit
# hanger positions are selected); re-assigning keeps every input across mode switches.
PERSISTENT = ["name", "description", "branch_participation", "m_length", "m_npipes", "b_length", "b_npipes",
              *(f"{p}_{k}" for p in ("m", "p", "n") for k in ("hmode", "hfirst", "hspacing", "hclear", "hpos")),
              "p_tees", "n_tees",
              "s_explicit", "s_values", "s_start", "s_stop", "s_n", "s_warm", "s_maxit", "s_tol", "s_test", "s_stol",
              "s_smaxit", "s_sdof_dc", "s_sdof_grid", "s_split", "th_xi", "th_test", "th_tol", "th_maxit", "v3_xi",
              "v3_test", "v3_tol", "v3_maxit", "v3_fb_test", "v3_fb_tol", "m_set", "m_levels", "m_floor",
              *(f"pipe_{k}" for k in PIPE_KEYS if k != "fluid_density")]
for _k in PERSISTENT:
    if _k in st.session_state:
        st.session_state[_k] = st.session_state[_k]


def parse_list(text: str) -> list[float]:
    return [float(v) for v in text.replace(";", ",").replace("\n", ",").split(",") if v.strip()]


# ============================================================================ sidebar
with st.sidebar:
    st.header("Inputs")
    src = st.radio("Start from", ["Archetype", "Example", "Upload file"], horizontal=True)
    if src == "Upload file":
        up = st.file_uploader("Layout file (YAML)", type=["yaml", "yml"])
        st.caption("Trapeze files referenced by the YAML can't be resolved from an upload; upload custom "
                   "trapezes in the Trapezes tab instead.")
        if up is not None and st.button("Load layout", type="primary"):
            try:
                data = yaml.safe_load(up.getvalue()) or {}
                for kind in KINDS:   # referenced files are not available
                    if isinstance((data.get("trapezes") or {}).get(kind), str):
                        data["trapezes"][kind] = "default"
                load_into_state(layout_from_dict(data))
                st.rerun()
            except Exception as exc:  # noqa: BLE001
                st.error(f"Could not load: {exc}")
    else:
        files = sorted(INPUT_DIRS[src].glob("*.yaml"))
        choice = st.selectbox(src, [f.stem for f in files])
        if src == "Archetype":
            st.caption("The paper's archetypes, read from their 3D models. For M01x–M03x, M61x–M63x and M61y the "
                       "paper's 2D models differ from these (docs/legacy_issues.md, E6).")
        if st.button("Load", type="primary"):
            try:
                load_into_state(load_layout_file(INPUT_DIRS[src] / f"{choice}.yaml"))
                st.rerun()
            except Exception as exc:  # noqa: BLE001
                st.error(f"Could not load: {exc}")

    st.divider()
    up_s = st.file_uploader("Settings file (YAML, optional)", type=["yaml", "yml"])
    c1, c2 = st.columns(2)
    if c1.button("Load settings", disabled=up_s is None):
        try:
            load_settings_into_state(AnalysisSettings.from_dict(yaml.safe_load(up_s.getvalue()) or {}))
            st.rerun()
        except Exception as exc:  # noqa: BLE001
            st.error(f"Could not load settings: {exc}")
    if c2.button("Default settings"):
        load_settings_into_state(AnalysisSettings.from_dict())
        st.rerun()
    st.divider()
    st.caption("Units: N, mm, s. Masses in tonnes.")

ss = st.session_state
v = ss.v
st.title("Suspended piping systems: reduced-order model")
tab_geo, tab_trap, tab_set, tab_po, tab_sdof, tab_th, tab_3d = st.tabs(
    ["Layout", "Trapezes", "Settings", "Pseudo-pushover", "SDOF parameters", "SDOF time history",
     "3D verification"])


# ============================================================================ layout inputs
def hanger_spec(prefix: str) -> dict:
    if ss[f"{prefix}_hmode"] == "Explicit positions":
        try:
            return {"positions": parse_list(ss[f"{prefix}_hpos"])}
        except ValueError:
            raise InputError("Hanger positions must be numbers separated by commas") from None
    return {"first": ss[f"{prefix}_hfirst"], "spacing": ss[f"{prefix}_hspacing"],
            "end_clearance": ss[f"{prefix}_hclear"]}


def hanger_positions(prefix: str, length: float) -> np.ndarray:
    """Hangers of a line with the layout's rules (empty if the inputs are invalid)."""
    from piperom.layout import _hanger_positions
    return _hanger_positions("hangers", hanger_spec(prefix), length, 1000.0)


def tees(side: str) -> list[float]:
    try:
        return parse_list(ss[f"{SIDE_KEY[side]}_tees"])
    except ValueError:
        raise InputError(f"Tee positions of the {side} branches must be numbers separated by commas") from None


def layout_dict() -> dict:
    """The layout description of the current inputs (piperom.layout format)."""
    pipe = {k: ss[f"pipe_{k}"] for k in PIPE_KEYS}
    main = {"length": ss.m_length, "n_pipes": int(ss.m_npipes), "hangers": hanger_spec("m"),
            "braces": {"transverse": {"positions": sorted(ss.m_brace_t)},
                       "longitudinal": {"positions": sorted(ss.m_brace_l)}}}
    branches = []
    for side, s in SIDE_KEY.items():
        for x in tees(side):
            branches.append({"x": x, "side": side, "length": ss.b_length, "n_pipes": int(ss.b_npipes),
                             "hangers": hanger_spec(s),
                             "braces": {"transverse": {"positions": sorted(ss[f"{s}_brace_t"])},
                                        "longitudinal": {"positions": sorted(ss[f"{s}_brace_l"])}}})
    return {"name": ss.name, "description": ss.description, "pipe": pipe, "main_line": main, "branches": branches}


def layout_from_state() -> tuple[Layout | None, dict[str, PipingSystem] | None, str | None]:
    try:
        lay = layout_from_dict(layout_dict())
        systems = lay.systems()
        for s in systems.values():
            s.trapezes = dict(ss.get("trapezes", {}))
            s.branch_participation = ss.branch_participation
            s.resolve()
        return lay, systems, None
    except (InputError, ValueError) as exc:
        return None, None, str(exc)


def hanger_inputs(prefix: str, length: float, title: str) -> np.ndarray:
    st.radio(f"{title}: hangers", ["Regular grid", "Explicit positions"], key=f"{prefix}_hmode", horizontal=True)
    if ss[f"{prefix}_hmode"] == "Regular grid":
        g1, g2, g3 = st.columns(3)
        g1.number_input("First (mm)", min_value=0.0, step=100.0, key=f"{prefix}_hfirst")
        g2.number_input("Spacing (mm)", min_value=1.0, step=100.0, key=f"{prefix}_hspacing")
        g3.number_input("End clearance (mm)", min_value=0.0, step=100.0, key=f"{prefix}_hclear")
    else:
        st.text_area("Positions (mm, comma separated)", key=f"{prefix}_hpos", height=68)
    try:
        return hanger_positions(prefix, length)
    except (InputError, ValueError) as exc:
        st.error(str(exc))
        return np.array([])


def brace_table(prefix: str, hx: np.ndarray, across: str, along: str) -> None:
    """One row per hanger with the two restraint directions as checkboxes."""
    bt, bl = f"{prefix}_brace_t", f"{prefix}_brace_l"
    t, l_ = {round(x, 6) for x in ss[bt]}, {round(x, 6) for x in ss[bl]}
    df = pd.DataFrame({"hanger": np.arange(1, len(hx) + 1), "position (mm)": hx,
                       across: [round(x, 6) in t for x in hx], along: [round(x, 6) in l_ for x in hx]})
    key = hashlib.md5(hx.tobytes()).hexdigest()[:8]
    edited = st.data_editor(df, key=f"braces_{prefix}_{v}_{key}", hide_index=True, height=min(250, 38 + 35 * len(hx)),
                            disabled=["hanger", "position (mm)"], width="stretch")
    if len(hx):   # keep the restraints while the hanger list is invalid
        ss[bt] = {float(x) for x, b in zip(edited["position (mm)"], edited[across]) if b}
        ss[bl] = {float(x) for x, b in zip(edited["position (mm)"], edited[along]) if b}


def plan_view(lay: Layout) -> go.Figure:
    """The layout in plan: lines, gravity hangers and restraints in each direction."""
    fig = go.Figure()
    L = lay.main.length / 1e3

    def add(xs, ys, name, marker, group, show):
        if len(xs):
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", marker=marker, name=name, legendgroup=group,
                                     showlegend=show))

    shown = set()

    def add_line_points(x, y, mask_y, mask_x):
        for sel, name, marker in [
                ((mask_y == 0) & (mask_x == 0), "Gravity hanger", dict(symbol="circle-open", size=7, color=C_MUTED)),
                ((mask_y == 1) & (mask_x == 0), "Restrained in y", dict(symbol="triangle-up", size=11, color=C_BRACE)),
                ((mask_y == 0) & (mask_x == 1), "Restrained in x", dict(symbol="triangle-right", size=11, color=C_BRACE_X)),
                ((mask_y == 1) & (mask_x == 1), "Restrained in x and y", dict(symbol="star", size=13, color="#ff7f0e"))]:
            add(x[sel], y[sel], name, marker, name, bool(name not in shown and sel.any()))
            if sel.any():
                shown.add(name)

    fig.add_trace(go.Scatter(x=[0, L], y=[0, 0], mode="lines", line=dict(color=C_MAIN, width=4), name="Main line"))
    m = lay.main
    add_line_points(m.hangers / 1e3, 0 * m.hangers, m.transverse, m.longitudinal)
    for k, b in enumerate(lay.branches):
        sgn = 1 if b.side == "+y" else -1
        x = b.x / 1e3
        fig.add_trace(go.Scatter(x=[x, x], y=[0, sgn * b.line.length / 1e3], mode="lines",
                                 line=dict(color=C_BRANCH, width=3), name="Branch", legendgroup="branch",
                                 showlegend=k == 0, hoverinfo="text",
                                 hovertext=f"Branch at x = {x:g} m ({b.side}), {b.line.n_pipes} pipes"))
        h = b.line.hangers / 1e3
        # along a branch (y): transverse restraints act along x, longitudinal along y
        add_line_points(np.full(len(h), x), sgn * h, b.line.longitudinal, b.line.transverse)
    fig.update_layout(height=520, margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.08),
                      xaxis_title="x (m)", yaxis_title="y (m)", title="Plan view")
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    return fig


def equivalent_view(system: PipingSystem, direction: str) -> go.Figure:
    """The equivalent line of one loading direction: hangers, restraints and branch DOFs."""
    rs = system.resolve()
    fig = go.Figure()
    L = rs.length / 1e3
    fig.add_trace(go.Scatter(x=[0, L], y=[0, 0], mode="lines", line=dict(color=C_MAIN, width=3), showlegend=False,
                             hoverinfo="skip"))
    hx = rs.hanger_x / 1e3
    fig.add_trace(go.Scatter(x=hx[rs.brace_mask == 0], y=0 * hx[rs.brace_mask == 0], mode="markers", showlegend=False,
                             marker=dict(symbol="circle-open", size=7, color=C_MUTED), hoverinfo="x"))
    fig.add_trace(go.Scatter(x=hx[rs.brace_mask == 1], y=0 * hx[rs.brace_mask == 1], mode="markers", showlegend=False,
                             marker=dict(symbol="triangle-up", size=11, color=C_BRACE), hoverinfo="x"))
    for j in range(rs.n_branches):
        x = rs.branch_x[j] / 1e3
        fig.add_trace(go.Scatter(x=[x, x], y=[0, 0.35], mode="lines+markers", showlegend=False,
                                 line=dict(color=C_BRANCH, width=2), marker=dict(symbol="diamond", size=[0, 10]),
                                 hoverinfo="text", hovertext=f"x = {x:g} m: {rs.branch_n_pipes[j]} pipes, "
                                                            f"L = {rs.branch_length[j] / 1e3:g} m, "
                                                            f"{rs.branch_n_braces[j]} braces"))
        fig.add_annotation(x=x, y=0.35, text=f"{rs.branch_n_braces[j]}", showarrow=False, yshift=10,
                           font=dict(color=C_BRANCH, size=10))
    if rs.n_mains_left == rs.n_mains_right:
        if rs.n_mains_left > 1:
            fig.add_annotation(x=L / 2, y=-0.3, text=f"{rs.n_mains_left} lines lumped", showarrow=False)
    else:
        xc = rs.x_center / 1e3
        fig.add_vline(x=xc, line=dict(color=C_MUTED, dash="dot"))
        fig.add_annotation(x=xc / 2, y=-0.3, text=f"{rs.n_mains_left} lines", showarrow=False)
        fig.add_annotation(x=(xc + L) / 2, y=-0.3, text=f"{rs.n_mains_right} lines", showarrow=False)
    what = ("a branch line (the branches lumped), the main line as its branch DOF" if direction == "x"
            else "the main line, the branches as DOFs")
    fig.update_layout(height=240, margin=dict(l=10, r=10, t=50, b=10), title=f"Loading along {direction}: {what}",
                      title_font=dict(size=12), xaxis_title="position along the analysed line (m)",
                      yaxis=dict(visible=False, range=[-0.6, 0.7]))
    return fig


with tab_geo:
    c_in, c_plot = st.columns([1, 2])
    with c_in:
        st.text_input("Name", key="name")
        st.text_area("Description", key="description", height=68)
        with st.expander("Pipe properties"):
            st.number_input("Outer diameter (mm)", min_value=0.1, key="pipe_outer_diameter")
            st.number_input("Inner diameter (mm)", min_value=0.0, key="pipe_inner_diameter")
            st.number_input("Elastic modulus (MPa)", min_value=1.0, key="pipe_elastic_modulus")
            st.number_input("Shear modulus (MPa)", min_value=1.0, key="pipe_shear_modulus")
            st.number_input("Density (t/mm³)", min_value=0.0, format="%.3e", key="pipe_density")
            fluid_default = st.checkbox("Fluid density = density / 7.8 (water in steel)",
                                        value=ss.pipe_fluid_density is None, key=f"fluid_default_{v}")
            if fluid_default:
                ss.pipe_fluid_density = None
            else:
                ss.pipe_fluid_density = st.number_input(
                    "Fluid density (t/mm³)", min_value=0.0, format="%.3e",
                    value=float(ss.pipe_fluid_density if ss.pipe_fluid_density is not None
                                else ss.pipe_density / 7.8), key=f"fluid_{v}")
            st.number_input("Mass factor", min_value=0.01, key="pipe_mass_factor")
            st.number_input("Branch participation factor", min_value=0.01, key="branch_participation")

        st.subheader("Main line (along x)")
        m1, m2 = st.columns(2)
        m1.number_input("Length (mm)", min_value=1.0, step=500.0, key="m_length")
        m2.number_input("Number of pipes", min_value=1, step=1, key="m_npipes")
        hx_main = hanger_inputs("m", ss.m_length, "Main line")
        st.caption("Restraints of the main line: transverse = along y, longitudinal = along x.")
        brace_table("m", hx_main, "restrained in y", "restrained in x")

        st.subheader("Branches (along y)")
        st.caption("All branches have the same length and pipes; the branches on each side of the main line are "
                   "identical (hangers and restraints, measured from the tee). Tees can be anywhere on the main line.")
        b1, b2 = st.columns(2)
        b1.number_input("Length (mm)", min_value=1.0, step=500.0, key="b_length")
        b2.number_input("Number of pipes", min_value=1, step=1, key="b_npipes")
        for side, s in SIDE_KEY.items():
            with st.expander(f"Branches on the {side} side", expanded=bool(ss[f"{s}_tees"].strip())):
                st.text_input("Tee positions along the main line (mm, comma separated; empty = none)",
                              key=f"{s}_tees")
                if ss[f"{s}_tees"].strip():
                    hx_b = hanger_inputs(s, ss.b_length, f"{side} branches")
                    brace_table(s, hx_b, "restrained in x", "restrained in y")



# ============================================================================ trapezes
def trapeze_figure(p: Pinching4, title: str) -> go.Figure:
    tri = trilinear_from_pinching4(p)
    d = [0.0, *p.pos_disp]
    f = [0.0, *p.pos_force]
    dn, fn = [0.0, *p.neg_disp], [0.0, *p.neg_force]
    u = np.linspace(0, p.pos_disp[-1], 200)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=d, y=np.array(f) / 1e3, mode="lines+markers", name="Pinching4 envelope",
                             line=dict(color=C_MAIN)))
    fig.add_trace(go.Scatter(x=dn, y=np.array(fn) / 1e3, mode="lines+markers", name="(negative)",
                             line=dict(color=C_MAIN, dash="dot"), showlegend=False))
    fig.add_trace(go.Scatter(x=u, y=[tri.force(x) / 1e3 for x in u], mode="lines",
                             name="Trilinear used in static procedure", line=dict(color=C_BRACE, dash="dash")))
    fig.update_layout(title=title, height=330, margin=dict(l=10, r=10, t=40, b=10),
                      xaxis_title="Deformation (mm)", yaxis_title="Force (kN)", legend=dict(orientation="h", y=-0.25))
    return fig


with tab_trap:
    st.caption("All braces of a type share one definition. Defaults: `inputs/trapezes/`. "
               "Upload a file in the same format or edit the values directly.")
    ss.trapezes = {}
    cols = st.columns(2)
    for col, (kind, label) in zip(cols, KINDS.items()):
        with col:
            st.subheader(label)
            st.caption(f"Source: {ss.trapeze_source[kind]}")
            up = st.file_uploader(f"Upload {kind} trapeze (CSV)", type=["csv"], key=f"up_{kind}_{v}")
            b1, b2 = st.columns(2)
            if b1.button("Use uploaded file", key=f"use_{kind}", disabled=up is None):
                try:
                    ss.trapeze_base[kind] = parse_trapeze_csv(up.getvalue().decode(), name=up.name)
                    ss.trapeze_source[kind] = up.name
                    ss.v += 1
                    st.rerun()
                except ValueError as exc:
                    st.error(str(exc))
            if b2.button("Reset to default", key=f"reset_{kind}"):
                ss.trapeze_base[kind] = load_trapeze("default", kind)
                ss.trapeze_source[kind] = "default file"
                ss.v += 1
                st.rerun()
            base = ss.trapeze_base[kind]
            params = base.parameters()
            df = pd.DataFrame({"parameter": PARAMETER_NAMES, "value": [str(params[k]) for k in PARAMETER_NAMES]})
            with st.expander("Parameters"):
                edited = st.data_editor(df, key=f"trap_{kind}_{v}", hide_index=True, disabled=["parameter"],
                                        width="stretch", height=400)
            try:
                p = Pinching4.from_parameters(dict(zip(edited["parameter"], edited["value"])), name=base.name)
                if p.parameters() != base.parameters():
                    ss.trapeze_source[kind] = "edited in app"
                ss.trapezes[kind] = p
                st.plotly_chart(trapeze_figure(p, label), width="stretch")
                st.download_button("Download CSV", p.to_csv(PREFIX[kind]), file_name=f"trapeze_{kind}.csv",
                                   key=f"dl_{kind}")
            except ValueError as exc:
                st.error(str(exc))

layout, systems, error = layout_from_state()     # again, with the trapezes of this run


# ============================================================================ settings
def settings_from_state() -> tuple[AnalysisSettings | None, str | None]:
    try:
        values = parse_list(ss.s_values) if ss.s_explicit else None
        data = {"pushover": {"delta_c_start": ss.s_start, "delta_c_stop": ss.s_stop, "n_steps": ss.s_n,
                             "delta_c_values": values, "warm_start": ss.s_warm},
                "shape_iteration": {"max_iterations": ss.s_maxit, "tolerance": ss.s_tol},
                "static_solver": {"test": ss.s_test, "tolerance": ss.s_stol, "max_iterations": ss.s_smaxit},
                "equivalent_static": {"branch_split": ss.s_split},
                "sdof": {**ss.s_raw["sdof"], "delta_c": ss.s_sdof_dc, "on_pushover_grid": ss.s_sdof_grid},
                "motions": {**ss.s_raw["motions"], "set": ss.m_set, "levels": list(ss.m_levels) or None,
                            "floor": int(ss.m_floor)},
                "sdof_time_history": {**ss.s_raw["sdof_time_history"], "damping_ratio": ss.th_xi,
                                      "test": ss.th_test, "tolerance": ss.th_tol, "max_iterations": ss.th_maxit},
                "verification_3d": {**ss.s_raw["verification_3d"], "damping_ratio": ss.v3_xi, "test": ss.v3_test,
                                    "tolerance": ss.v3_tol, "max_iterations": ss.v3_maxit,
                                    "fallback_test": ss.v3_fb_test, "fallback_tolerance": ss.v3_fb_tol}}
        return AnalysisSettings.from_dict(data), None
    except (InputError, ValueError) as exc:
        return None, str(exc)


with tab_set:
    c1, c2, c3 = st.columns(3)
    with c1:
        st.subheader("Pseudo-pushover")
        st.checkbox("Explicit list of target displacements", key="s_explicit")
        if ss.s_explicit:
            st.text_area("Delta_c values (mm)", key="s_values")
        else:
            st.number_input("Delta_c start (mm)", min_value=1e-6, format="%.3f", key="s_start")
            st.number_input("Delta_c stop (mm)", min_value=1e-6, format="%.3f", key="s_stop")
            st.number_input("Number of steps", min_value=1, step=1, key="s_n")
        st.checkbox("Warm start from previous step's shape", key="s_warm")
    with c2:
        st.subheader("Shape iteration")
        st.number_input("Maximum iterations", min_value=1, step=1, key="s_maxit")
        st.number_input("Tolerance (RMS shape change)", min_value=1e-12, format="%.1e", key="s_tol")
        st.subheader("Static solver")
        st.selectbox("Convergence test", ["NormDispIncr", "EnergyIncr", "NormUnbalance"], key="s_test")
        st.number_input("Solver tolerance", min_value=1e-16, format="%.1e", key="s_stol")
        st.number_input("Solver max iterations", min_value=1, step=1, key="s_smaxit")
        st.selectbox("Branch force split", ["consistent", "legacy"], key="s_split",
                     help="consistent: the branch force is fully applied (junction node mass in the main-line "
                          "share), as the paper's scripts. legacy: as the original code before the A2 fix "
                          "(docs/legacy_issues.md).")
    with c3:
        st.subheader("Equivalent SDOF")
        st.number_input("Delta_c defining the SDOF (mm)", min_value=1e-6, format="%.3f", key="s_sdof_dc")
        st.checkbox("Use the closest pushover step", key="s_sdof_grid",
                    help="As Pushover_SDOF/sdof_from_2d.py, which reads the SDOF from the pushover results.")
    tests = ["EnergyIncr", "RelativeEnergyIncr", "NormDispIncr", "NormUnbalance"]
    c4, c5 = st.columns(2)
    with c4:
        st.subheader("SDOF time history")
        st.number_input("Damping ratio", min_value=0.0, max_value=1.0, format="%.3f", key="th_xi")
        st.selectbox("Convergence test", tests, key="th_test")
        st.number_input("Tolerance", min_value=1e-16, format="%.1e", key="th_tol")
        st.number_input("Max iterations", min_value=1, step=1, key="th_maxit")
        st.caption("Mass-proportional damping at the initial SDOF frequency, springs capped at the 3D restraint "
                   "caps, stop when a braced node exceeds the collapse displacement (as Pushover_SDOF/*_NLTHA.py).")
    with c5:
        st.subheader("3D verification")
        st.number_input("Damping ratio ", min_value=0.0, max_value=1.0, format="%.3f", key="v3_xi")
        st.selectbox("Convergence test ", tests, key="v3_test")
        st.number_input("Tolerance ", min_value=1e-16, format="%.1e", key="v3_tol")
        st.number_input("Max iterations ", min_value=1, step=1, key="v3_maxit")
        st.selectbox("Fallback test", ["script"] + tests, key="v3_fb_test",
                     help="script: the test types of the model's own 3D script (they differ between models).")
        st.number_input("Fallback tolerance", min_value=1e-16, format="%.1e", key="v3_fb_tol")
        st.caption("Mass-proportional damping at mode 1; stop when a braced node exceeds the collapse "
                   "displacement (as the paper's 3D_models/*_biron.py).")
    settings, s_error = settings_from_state()
    if s_error:
        st.error(s_error)
    else:
        st.download_button("Download settings (YAML)", settings.to_yaml(), file_name="settings.yaml")


# ============================================================================ bundles
def custom_trapeze_files() -> dict[str, str]:
    out = {}
    for kind, p in ss.get("trapezes", {}).items():
        if isinstance(p, Pinching4) and p.parameters() != load_trapeze("default", kind).parameters():
            out[f"custom_{kind}.csv"] = p.to_csv(PREFIX[kind])
    return out


def trapeze_refs() -> dict[str, str]:
    custom = custom_trapeze_files()
    return {k: (f"custom_{k}.csv" if f"custom_{k}.csv" in custom else "default") for k in KINDS}


def layout_yaml() -> str:
    d = layout_dict()
    d["trapezes"] = trapeze_refs()
    return dump_yaml(d)


def system_yaml(system: PipingSystem) -> str:
    d = system.to_dict()
    d["trapezes"] = trapeze_refs()
    return dump_yaml(d)


def zip_bytes(files: dict[str, str]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in files.items():
            z.writestr(name, text)
    return buf.getvalue()


def input_files(settings: AnalysisSettings) -> dict[str, str]:
    return {"layout.yaml": layout_yaml(), **{f"system_{d}.yaml": system_yaml(systems[d]) for d in DIRS},
            "settings.yaml": settings.to_yaml(), **custom_trapeze_files()}


def signature(*parts: str) -> str:
    return hashlib.sha1("\n".join(parts).encode()).hexdigest()


with tab_geo:
    with c_plot:
        if error:
            st.error(error)
        else:
            p1, p2 = st.columns([3, 2])
            p1.plotly_chart(plan_view(layout), width="stretch")
            with p2:
                for d in DIRS:
                    st.plotly_chart(equivalent_view(systems[d], d), width="stretch", key=f"eq_{d}")
                st.caption("Equivalent lines analysed by the reduced-order model: ▲ restraint across the line, "
                           "○ gravity hanger, ◆ branch DOF (number of longitudinal restraints above it).")
            if settings is not None:
                st.download_button("Download inputs (layout + equivalent systems + settings + custom trapezes, zip)",
                                   zip_bytes(input_files(settings)), file_name=f"{ss.name}_inputs.zip")

ready = systems is not None and settings is not None
po_sig, sd_sig = {}, {}
if ready:
    trap_sig = "".join(str(sorted(p.parameters().items())) for p in ss.trapezes.values())
    base_set = dump_yaml({k: val for k, val in settings.to_dict().items()
                          if k not in ("sdof", "motions", "sdof_time_history", "verification_3d")})
    for d in DIRS:
        po_sig[d] = signature(system_yaml(systems[d]), trap_sig, base_set)
        sd_sig[d] = po_sig[d] + signature(dump_yaml(settings.to_dict()["sdof"]))


def not_ready_warning() -> None:
    st.warning("Fix the inputs first: " + (error or s_error or ""))


# ============================================================================ SDOF backbone
def sdof_backbone(params) -> tuple[np.ndarray, np.ndarray]:
    """Monotonic backbone of the SDOF: sum of the springs' positive envelopes, up to the first envelope end."""
    u_max = min(s.spring.pos_disp[-1] for s in params.supports)
    u = np.linspace(0, u_max, 300)
    F = sum(np.interp(u, [0, *s.spring.pos_disp], [0, *s.spring.pos_force]) for s in params.supports)
    return u, F


# ============================================================================ pushover
def run_remote(fn, *args, label="Step"):
    """Run ``fn(*args, queue)`` in a worker process, showing its (k, n) progress messages."""
    q = manager().Queue()
    fut = executor().submit(fn, *args, q)
    bar = st.progress(0.0, text="Starting worker...")
    while not fut.done():
        try:
            k, n = q.get(timeout=0.2)
            text = f"{label} {k}/{n}" if isinstance(k, int) else f"{label} {100 * k / n:.0f}%"
            bar.progress(min(k / n, 1.0), text=text)
        except Exception:  # noqa: BLE001  (queue.Empty)
            pass
    bar.empty()
    return fut.result()


def show_pushover(d: str) -> None:
    res = ss.po_result.get(d)
    if res is None:
        st.info("Not run yet.")
        return
    if ss.po_sig.get(d) != po_sig[d]:
        st.info("Inputs changed since this run; results below refer to the previous inputs.")
    if res.n_not_converged:
        st.warning(f"{res.n_not_converged} step(s) did not reach the shape tolerance.")
    curve = pd.DataFrame(res.curve_rows())
    c1, c2 = st.columns(2)
    with c1:
        xaxis = st.radio("Capacity curve abscissa", ["u_sdof", "delta_c"], horizontal=True, key=f"po_x_{d}",
                         format_func={"u_sdof": "Equivalent SDOF displacement", "delta_c": "Target displacement Δc"}.get)
        fig = go.Figure(go.Scatter(x=curve[xaxis], y=curve.base_shear / 1e3, mode="lines+markers",
                                   line=dict(color=C_MAIN), name="2D pseudo-pushover"))
        sd = ss.sdof_result.get(d)
        if sd is not None and ss.sdof_sig.get(d, "").startswith(po_sig[d]) and xaxis == "u_sdof":
            u, F = sdof_backbone(sd)
            fig.add_trace(go.Scatter(x=u, y=F / 1e3, mode="lines", line=dict(color=C_BRACE, dash="dash"),
                                     name=f"SDOF backbone (Δc = {sd.delta_c:g} mm)"))
        fig.update_layout(height=380, margin=dict(l=10, r=10, t=10, b=10), yaxis_title="Base shear (kN)",
                          xaxis_title="u_SDOF (mm)" if xaxis == "u_sdof" else "Δc (mm)",
                          legend=dict(orientation="h", y=1.1))
        st.plotly_chart(fig, width="stretch", key=f"po_curve_{d}")
    with c2:
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=curve.delta_c, y=curve.gamma, name="Γ", line=dict(color=C_MAIN)))
        fig.add_trace(go.Scatter(x=curve.delta_c, y=curve.mass_ratio, name="Effective mass ratio",
                                 line=dict(color=C_BRANCH)))
        fig.update_layout(height=410, margin=dict(l=10, r=10, t=40, b=10), xaxis_title="Δc (mm)",
                          legend=dict(orientation="h", y=1.1))
        st.plotly_chart(fig, width="stretch", key=f"po_gamma_{d}")

    k = st.select_slider("Displaced shape at Δc (mm)", options=list(range(len(res.steps))), key=f"po_k_{d}",
                         format_func=lambda i: f"{res.steps[i].delta_c:.3g}", value=len(res.steps) - 1)
    step = res.steps[k]
    dofs = pd.DataFrame(res.dofs)
    scaled = st.toggle("Show displacements in mm (otherwise normalised to the reference branch DOF)",
                       key=f"po_scaled_{d}")
    yv = step.d_scaled if scaled else step.d_norm
    h = dofs.kind == "hanger"
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=dofs.x[h] / 1e3, y=yv[h.values], mode="lines+markers", name="Analysed line",
                             line=dict(color=C_MAIN), marker=dict(size=5)))
    hb = h & dofs.braced
    fig.add_trace(go.Scatter(x=dofs.x[hb] / 1e3, y=yv[hb.values], mode="markers", name="Restraint",
                             marker=dict(symbol="triangle-up", size=11, color=C_BRACE)))
    b = dofs.kind == "branch"
    fig.add_trace(go.Scatter(x=dofs.x[b] / 1e3, y=yv[b.values], mode="markers", name="Branch DOF",
                             marker=dict(symbol="diamond", size=11, color=C_BRANCH),
                             hovertext=[f"branch {int(j) + 1}" for j in dofs.branch[b]]))
    fig.update_layout(height=360, margin=dict(l=10, r=10, t=10, b=10), xaxis_title="position along the analysed line (m)",
                      yaxis_title="Displacement (mm)" if scaled else "Normalised displacement",
                      legend=dict(orientation="h", y=1.1))
    st.plotly_chart(fig, width="stretch", key=f"po_shape_{d}")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Base shear", f"{step.base_shear / 1e3:.2f} kN")
    m2.metric("Γ", f"{step.gamma:.4f}")
    m3.metric("Effective mass", f"{step.effective_mass:.3f} t")
    m4.metric("Shape iterations", f"{step.iterations}" + ("" if step.converged else " (not converged)"))
    with st.expander("Capacity curve table"):
        st.dataframe(curve, hide_index=True, width="stretch")
    files = {"pushover_curve.csv": rows_to_csv(res.curve_rows(), CURVE_COLUMNS),
             "pushover_shapes.csv": rows_to_csv(res.shape_rows(), SHAPE_COLUMNS), **ss.po_inputs}
    d1, d2 = st.columns(2)
    d1.download_button("Download capacity curve (CSV)", files["pushover_curve.csv"],
                       file_name=f"{res.name}_pushover_curve.csv", key=f"po_dl1_{d}")
    d2.download_button("Download all results and inputs (zip)", zip_bytes(files),
                       file_name=f"{res.name}_pushover.zip", key=f"po_dl2_{d}")


for _key in ("po_result", "po_sig", "sdof_result", "sdof_sig", "th_result", "th_sig"):
    ss.setdefault(_key, {})

with tab_po:
    if not ready:
        not_ready_warning()
    else:
        n_steps = len(settings.delta_c())
        st.write(f"**{ss.name}**, both loading directions: {n_steps} target displacements from "
                 f"{settings.delta_c()[0]:g} to {settings.delta_c()[-1]:g} mm.")
        if st.button("Run pseudo-pushover", type="primary"):
            try:
                for d in DIRS:
                    ss.po_result[d] = run_remote(pushover_job, systems[d], settings, label=f"Loading along {d}: step")
                    ss.po_sig[d] = po_sig[d]
                ss.po_inputs = input_files(settings)
            except Exception as exc:  # noqa: BLE001
                st.error(f"Analysis failed: {exc}")
        for d, t in zip(DIRS, st.tabs([f"Loading along {d}" for d in DIRS])):
            with t:
                show_pushover(d)


# ============================================================================ SDOF
def current_sdof(d: str):
    """SDOF parameters derived in this session for the current inputs, or None."""
    sd = ss.sdof_result.get(d)
    return sd if ready and sd is not None and ss.sdof_sig.get(d) == sd_sig[d] else None


def derive_sdof_button(key: str) -> None:
    if st.button("Derive SDOF parameters", type="primary", key=key):
        try:
            with st.spinner("Running the equivalent static procedure (both directions)..."):
                for d in DIRS:
                    ss.sdof_result[d] = executor().submit(sdof_job, systems[d], settings, settings.sdof_delta_c).result()
                    ss.sdof_sig[d] = sd_sig[d]
            ss.sdof_inputs = input_files(settings)
        except Exception as exc:  # noqa: BLE001
            st.error(f"Derivation failed: {exc}")
        else:
            st.rerun()


def show_sdof(d: str) -> None:
    sd = ss.sdof_result.get(d)
    if sd is None:
        st.info("Not derived yet.")
        return
    if ss.sdof_sig.get(d) != sd_sig[d]:
        st.info("Inputs changed since these parameters were derived.")
    if not sd.converged:
        st.warning("The displaced shape did not reach the tolerance at this Δc.")
    m = st.columns(5)
    m[0].metric("Γ", f"{sd.gamma:.4f}")
    m[1].metric("Effective mass", f"{sd.effective_mass:.3f} t")
    m[2].metric("Mass ratio", f"{sd.mass_ratio:.3f}")
    m[3].metric("u_SDOF at Δc", f"{sd.u_sdof:.3f} mm")
    m[4].metric("Base shear at Δc", f"{sd.base_shear / 1e3:.2f} kN")
    rows = [{"kind": s.kind, "position (mm)": s.x, "branch": "" if s.branch is None else str(s.branch + 1),
             "φ": s.phi, "trapezes": s.n_trapezes, "Γ·φ": sd.gamma * s.phi,
             **{f"ePd{i + 1} (mm)": s.spring.pos_disp[i] for i in range(4)},
             **{f"ePf{i + 1} (kN)": s.spring.pos_force[i] / 1e3 for i in range(4)}}
            for s in sd.supports]
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    u, F = sdof_backbone(sd)
    fig = go.Figure(go.Scatter(x=u, y=F / 1e3, mode="lines", line=dict(color=C_BRACE), name="SDOF backbone"))
    fig.update_layout(height=320, margin=dict(l=10, r=10, t=10, b=10), xaxis_title="u_SDOF (mm)",
                      yaxis_title="Force (kN)")
    st.plotly_chart(fig, width="stretch", key=f"sd_backbone_{d}")
    st.caption("Backbone = sum of the scaled spring envelopes (positive branch), plotted up to the first "
               "spring's last envelope point. Overlay it on the pushover curve in the Pseudo-pushover tab.")
    files = {"sdof_parameters.yaml": sd.to_yaml(), "sdof_springs.csv": sd.springs_csv(), **ss.sdof_inputs}
    d1, d2 = st.columns(2)
    d1.download_button("Download SDOF parameters (YAML)", files["sdof_parameters.yaml"],
                       file_name=f"{sd.name}_sdof_parameters.yaml", key=f"sd_dl1_{d}")
    d2.download_button("Download all (zip)", zip_bytes(files), file_name=f"{sd.name}_sdof.zip", key=f"sd_dl2_{d}")


with tab_sdof:
    if not ready:
        not_ready_warning()
    else:
        st.write(f"Equivalent SDOF of **{ss.name}** in both directions at Δc = **{sdof_delta_c(settings):g} mm** "
                 "(change it in the Settings tab).")
        derive_sdof_button("derive_sdof")
        for d, t in zip(DIRS, st.tabs([f"Loading along {d}" for d in DIRS])):
            with t:
                show_sdof(d)


# ============================================================================ SDOF time history
def show_time_history(d: str, th_sig: str) -> None:
    if not ss.th_result.get(d):
        st.info("Not run yet.")
        return
    model, resp = ss.th_result[d]
    if ss.th_sig.get(d) != th_sig:
        st.info("Inputs or selection changed since this run; results refer to the previous run.")
    phi = max(model.support_phi)
    df = pd.DataFrame([{"level": r.level, "record": r.record, "peak u_SDOF (mm)": r.peak_u,
                        "peak support displacement (mm)": model.gamma * phi * r.peak_u,
                        "completed": r.completed, "collapsed": r.collapsed} for r in resp])
    if df.collapsed.any():
        st.info(f"{int(df.collapsed.sum())} analysis/es stopped at collapse (a braced node above "
                f"{settings.sdof_time_history['collapse_displacement']:g} mm).")
    bad = int((~df.completed & ~df.collapsed).sum())
    if bad:
        st.warning(f"{bad} analysis/es stopped before the end of the record.")
    ycol = "peak support displacement (mm)"
    fig = go.Figure()
    if df.level.notna().any() and df.level.nunique() > 1:
        fig.add_trace(go.Scatter(x=df.level, y=df[ycol], mode="markers", name="records",
                                 marker=dict(color=C_MAIN, opacity=0.45, size=7), text=df.record))
        med = df.groupby("level")[ycol].median()
        fig.add_trace(go.Scatter(x=med.index, y=med.values, mode="lines+markers", name="median",
                                 line=dict(color=C_BRACE, width=3)))
        fig.update_layout(xaxis_title="Intensity level")
    else:
        fig.add_trace(go.Bar(x=df.record.astype(str), y=df[ycol], marker_color=C_MAIN))
        fig.update_layout(xaxis_title="Record", xaxis_type="category")
    fig.update_layout(height=360, margin=dict(l=10, r=10, t=10, b=10), yaxis_title=ycol,
                      legend=dict(orientation="h", y=1.1))
    st.plotly_chart(fig, width="stretch", key=f"th_peaks_{d}")
    st.caption("Peak support displacement = Γ · max φ · peak u_SDOF (largest support demand).")
    with st.expander("Results table"):
        st.dataframe(df, hide_index=True, width="stretch")
    st.download_button("Download peaks (CSV)", df.to_csv(index=False), file_name=f"{model.name}_sdof_peaks.csv",
                       key=f"th_dl_{d}")
    with_hist = [i for i, r in enumerate(resp) if len(r.u)]
    if with_hist:
        i = st.selectbox("Response history", with_hist, key=f"th_hist_{d}",
                         format_func=lambda i: f"record {resp[i].record}"
                         + ("" if resp[i].level is None else f", IM{resp[i].level}"))
        r = resp[i]
        h1, h2 = st.columns(2)
        f1 = go.Figure(go.Scatter(x=r.time, y=r.u, line=dict(color=C_MAIN, width=1)))
        f1.update_layout(height=300, margin=dict(l=10, r=10, t=30, b=10), title="SDOF displacement",
                         xaxis_title="Time (s)", yaxis_title="u (mm)")
        h1.plotly_chart(f1, width="stretch", key=f"th_u_{d}")
        f2 = go.Figure(go.Scatter(x=r.u, y=r.force / 1e3, line=dict(color=C_BRACE, width=1)))
        f2.update_layout(height=300, margin=dict(l=10, r=10, t=30, b=10), title="Hysteresis",
                         xaxis_title="u (mm)", yaxis_title="Force (kN)")
        h2.plotly_chart(f2, width="stretch", key=f"th_f_{d}")
        st.dataframe(pd.DataFrame(support_demands(model, r.peak_u)), hide_index=True, width="stretch")
        st.caption("Support demands: peak displacement Γ·φ·u and its ratio to the trapeze envelope "
                   "deformations ePd1-ePd4.")
    else:
        st.caption("Response histories are kept for batches of up to 100 analyses.")


with tab_th:
    sets = load_motion_sets()
    if not ready:
        not_ready_warning()
    elif not sets:
        st.warning("No floor-motion sets: add them to motions/motion_sets.yaml.")
    elif any(current_sdof(d) is None for d in DIRS):
        st.info(f"The SDOF time history uses the SDOFs of the current inputs (Δc = {sdof_delta_c(settings):g} mm).")
        derive_sdof_button("derive_th")
    else:
        sds = {d: current_sdof(d) for d in DIRS}
        st.write(" · ".join(f"Loading along {d}: Γ = {sds[d].gamma:.4f}, M_eff = {sds[d].effective_mass:.3f} t, "
                            f"{len(sds[d].supports)} springs" for d in DIRS)
                 + ". The two directions are analysed independently, each under the selected records.")
        if ss.m_set not in sets:
            ss.m_set = next(iter(sets))
        c1, c2, c3 = st.columns([1, 1, 2])
        with c1:
            st.selectbox("Motion set", list(sets), key="m_set", help="Defined in motions/motion_sets.yaml")
            ms = sets[ss.m_set]
            st.caption(ms.description)
            floors = sorted(ms.floors)
            if ss.m_floor not in floors:
                ss.m_floor = floors[-1]
            st.selectbox("Floor", floors, key="m_floor")
        with c2:
            if ms.levels is not None:
                ss.m_levels = [lv for lv in ss.m_levels if lv in ms.levels] or [ms.levels[-1]]
                st.multiselect("Intensity levels", ms.levels, key="m_levels")
        with c3:
            all_rec = st.toggle("All records", value=True, key="th_all")
            recs = None if all_rec else st.multiselect("Records", ms.records, default=ms.records[:2], key="th_recs")
        try:
            runs = select_runs(ms, list(ss.m_levels) if ms.levels else None, recs)
        except InputError as exc:
            runs = []
            st.error(str(exc))
        st.caption(f"{len(runs)} analyses per direction (about 0.3 s each).")
        th_sig = {d: sd_sig[d] + signature(ms.name, str(runs), str(ss.m_floor), dump_yaml(settings.sdof_time_history))
                  for d in DIRS}
        if st.button("Run SDOF time histories", type="primary", disabled=not runs):
            try:
                for d in DIRS:
                    model = SDOFModel.from_parameters(sds[d])
                    ss.th_result[d] = (model, run_remote(sdof_time_history_job, model, ms.name, runs, int(ss.m_floor),
                                                         settings, label=f"Loading along {d}: analysis"))
                    ss.th_sig[d] = th_sig[d]
            except Exception as exc:  # noqa: BLE001
                st.error(f"Analysis failed: {exc}")
        for d, t in zip(DIRS, st.tabs([f"Loading along {d}" for d in DIRS])):
            with t:
                show_time_history(d, th_sig[d])


# ============================================================================ 3D verification
def node_plan(model3d, nodes, peaks, title):
    xy = model3d.node_coords()
    fig = go.Figure(go.Scatter(
        x=[xy[n][0] / 1e3 for n in nodes], y=[xy[n][1] / 1e3 for n in nodes], mode="markers",
        marker=dict(size=11, color=peaks, colorscale="Viridis", showscale=True,
                    colorbar=dict(title="mm", thickness=12)),
        text=[f"node {n}: {p:.2f} mm" for n, p in zip(nodes, peaks)], hoverinfo="text"))
    fig.update_layout(height=340, margin=dict(l=10, r=10, t=30, b=10), title=title, xaxis_title="x (m)",
                      yaxis_title="y (m)")
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    return fig


with tab_3d:
    models = list_models3d()
    paired = {n: m for n, m in load_motion_sets().items() if m.record_pairs}
    if not ready:
        not_ready_warning()
    elif not models or not paired:
        st.warning("Needs 3D models in inputs/models3d/ and a motion set with record pairs.")
    else:
        st.caption("Full 3D model of a paper archetype under a bidirectional floor motion, compared with the "
                   "reduced-order model of each direction (the paper's systems in inputs/archetypes/, SDOF at the "
                   "Δc of the Settings tab, default trapezes).")
        default = ss.name if ss.name in models else None
        c1, c2, c3, c4 = st.columns(4)
        name3d = c1.selectbox("3D model", models, index=models.index(default) if default in models else 0)
        set3d = c2.selectbox("Motion set ", list(paired))
        ms3 = paired[set3d]
        level3d = c3.selectbox("Intensity level", ms3.levels, index=len(ms3.levels) - 1) if ms3.levels else None
        floor3d = c4.selectbox("Floor ", sorted(ms3.floors), index=len(ms3.floors) - 1)
        c5, c6 = st.columns(2)
        pairs = ms3.pairs()
        k = c5.selectbox("Ground motion", range(len(pairs)), format_func=lambda i: f"{i + 1}: {pairs[i][0]} / {pairs[i][1]}")
        flip = c6.radio("Orientation", ["x ← first component", "x ← second component"], horizontal=True) != \
            "x ← first component"
        rx, ry = pairs[k][::-1] if flip else pairs[k]
        st.caption(f"x ← {rx}, y ← {ry}. One run takes about 45 s for M01–M03 and several minutes for the "
                   "larger models.")
        if st.button("Run 3D verification", type="primary"):
            try:
                ss.v3_result = run_remote(verification_job, name3d, set3d, rx, ry, level3d, int(floor3d), settings,
                                          label="3D analysis")
            except Exception as exc:  # noqa: BLE001
                st.error(f"Analysis failed: {exc}")
        res3 = ss.get("v3_result")
        if res3 is not None:
            r = res3.response
            m3 = load_model3d(r.model)
            (st.success if r.completed else st.error)(
                f"{r.model}: x ← {r.record_x}, y ← {r.record_y}"
                + ("" if r.level is None else f", IM{r.level}")
                + (" — completed." if r.completed else " — collapse: a braced node exceeded the collapse displacement "
                   f"at {r.end_time:.3f} s." if r.collapsed else f" — stopped at {r.end_time:.3f} of {r.duration:.3f} s."))
            rows = []
            at_braced = braced_peaks(m3, r)
            for d, peaks in (("x", r.peak_x), ("y", r.peak_y)):
                at_braces = at_braced[d]
                row = {"direction": d, "3D max peak, all nodes (mm)": peaks.max(),
                       "3D max peak at braced nodes (mm)": at_braces}
                p = res3.rom.get(d)
                if p is not None:
                    row.update({"ROM system": p.system, "Γ": p.gamma, "max φ": p.phi_max, "SDOF peak u (mm)": p.peak_u,
                                "ROM Γ·φ·u (mm)": p.peak_support_displacement,
                                "ROM / 3D (braced nodes)": p.peak_support_displacement / at_braces})
                rows.append(row)
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
            st.caption("'3D max peak at braced nodes' is the measure stored in the paper's Results/ "
                       "(DispX, DispY): the largest peak displacement over the braced nodes of the 3D model.")
            st.caption(f"3D periods: {', '.join(f'{T:.4f}' for T in r.periods)} s.")
            p1, p2 = st.columns(2)
            p1.plotly_chart(node_plan(m3, r.nodes_x, r.peak_x, "Peak displacement in x"), width="stretch")
            p2.plotly_chart(node_plan(m3, r.nodes_y, r.peak_y, "Peak displacement in y"), width="stretch")
            d = st.radio("History of node", ["x", "y"], horizontal=True, key="v3_dir")
            nodes, u = (r.nodes_x, r.ux) if d == "x" else (r.nodes_y, r.uy)
            peaks = r.peak_x if d == "x" else r.peak_y
            n = st.selectbox("Node", nodes, index=int(np.argmax(peaks)))
            f = go.Figure(go.Scatter(x=r.time, y=u[:, nodes.index(n)], line=dict(color=C_MAIN, width=1)))
            f.update_layout(height=300, margin=dict(l=10, r=10, t=10, b=10), xaxis_title="Time (s)",
                            yaxis_title=f"u_{d} (mm)")
            st.plotly_chart(f, width="stretch")
            csv = pd.DataFrame([{"direction": "x", "node": a, "peak": b} for a, b in zip(r.nodes_x, r.peak_x)]
                               + [{"direction": "y", "node": a, "peak": b} for a, b in zip(r.nodes_y, r.peak_y)])
            st.download_button("Download 3D peaks (CSV)", csv.to_csv(index=False), file_name=f"{r.model}_3d_peaks.csv")
