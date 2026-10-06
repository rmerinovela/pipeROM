"""
Full 3D OpenSees model of a suspended piping system, from a description of the whole system in plan.

The model is the one of generate_3d_model.py (the 3D models of the paper, generalized): a straight main
line along x from x = 0 to its length and straight branches along y framing into it, each a bundle of
identical pipes on gravity hangers. Restrained hangers carry trapeze springs (Pinching4) in the
transverse and/or longitudinal direction of their line, the others soft springs. Pipe joints are
threaded-joint hinges every 6 m along each line and at each tee.

Inputs use the conventions of the piperom system files (inputs/defaults/system.yaml, docs/input_files.md):
hangers as a grid or positions, restraints as a mask, positions or count. A count is placed with the app's
rule (spacing L/n starting half a spacing from the start of the line, snapped to the nearest free hanger),
on every line and in both directions; branch positions are measured from the tee.

The result is a dictionary in the format of inputs/models3d/*.json, run by piperom.verification3d:
  commands     OpenSees commands up to the gravity load pattern
  node_groups  X = Y: pipe nodes recorded for the peak and displaced shape (main line sorted by x,
               then each branch sorted by y); Xt / Yt: restrained in x / y; Xp / Yp: the others.

Description (YAML), units N, mm, s:

    name: MySystem
    description: free text                     # optional
    pipe: {...}                                # optional, as in the system files (defaults: inputs/defaults)
    main_line:
      length: 36000
      n_pipes: 3
      hangers: {first: 1000, spacing: 3000}    # or {positions: [...]}; optional end_clearance
      braces:
        transverse: {count: 4}                 # restraints along y; {mask: [...]} or {positions: [...]}
        longitudinal: {count: 3}               # restraints along x
    branches:
      - x: 36000                               # tee position on the main line
        side: +y                               # +y or -y
        length: 18000
        n_pipes: 3                             # optional, default: main line's
        hangers: {first: 1000, spacing: 3000}  # optional, default: main line's grid; distances from the tee
        braces:
          transverse: {count: 2}               # restraints along x
          longitudinal: {count: 2}             # restraints along y
    trapezes: {transverse: default, longitudinal: default}   # optional, as in the system files
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from piperom import DEFAULTS_DIR                                  # noqa: E402
from piperom.inputs import (InputError, check_keys, deep_merge, evenly_spaced_brace_mask,   # noqa: E402
                            parse_number, read_yaml)
from piperom.static_model import _section                         # noqa: E402
from piperom.trapeze import load_trapeze                          # noqa: E402

G_ACC = 9805.0              # mm/s^2
Z_HANGER, Z_SPRING = 10.0, 5.0
JOINT_SPACING = 6000.0      # threaded joints along each line
JOINT_HALF = 60.0           # half length of a joint (two hinges 120 mm apart)
BRANCH_START = 100.0        # first joint of a branch, from the tee
CAP_LONG, CAP_TRANS = 60.0, 35.0    # MinMax limits of the trapeze springs: end of each envelope
                                    # (61 mm longitudinal, 36 mm transverse); swapped w.r.t. the
                                    # paper's 3D models (docs/legacy_issues.md, D1)


# ----------------------------------------------------------------------------------------- inputs
def load_description(path: str | Path) -> dict:
    path = Path(path)
    d = read_yaml(path)
    d.setdefault("_base_dir", str(path.parent))
    return d


def _hanger_positions(section: str, spec: dict, length: float, defaults: dict) -> np.ndarray:
    """Hanger positions along a line (distance from its start), with the piperom rules."""
    check_keys(section, spec, {"first", "spacing", "end_clearance", "positions"})
    if spec.get("positions") is not None:
        x = np.round(np.array(spec["positions"], float), 6)
        if len(x) == 0 or np.any(np.diff(x) <= 0) or x[0] <= 0 or x[-1] >= length:
            raise InputError(f"{section}: positions must be increasing, strictly inside the line")
        return x
    x0 = parse_number(section, "first", spec.get("first"))
    s = parse_number(section, "spacing", spec.get("spacing"))
    c = parse_number(section, "end_clearance", spec.get("end_clearance", defaults["end_clearance"]), positive=False)
    if x0 > length - c:
        raise InputError(f"{section}: first hanger must be at least {c:g} mm from the end")
    n_max = int((length - c - x0) // s)
    return np.round(np.array([x0 + i * s for i in range(n_max + 1) if x0 + i * s <= length - c + 1e-6]), 6)


def _count_mask(hx: np.ndarray, length: float, n: int) -> np.ndarray:
    """The app's rule (piperom, compute_stiff_mask of the paper's equivalent_static.py): n restraints
    at length/n spacing starting half a spacing from the line start, each snapped to the nearest free
    hanger. One implementation for the app and the workflow."""
    return evenly_spaced_brace_mask(hx, length, n)


def _brace_mask(section: str, spec: dict | None, hx: np.ndarray, length: float) -> np.ndarray:
    spec = spec or {"count": 0}
    check_keys(section, spec, {"mask", "positions", "count"})
    if sum(k in spec for k in ("mask", "positions", "count")) != 1:
        raise InputError(f"{section}: give exactly one of mask, positions, count")
    if "mask" in spec:
        mask = np.array(spec["mask"], int)
        if len(mask) != len(hx) or not np.all((mask == 0) | (mask == 1)):
            raise InputError(f"{section}: mask needs one 0/1 per hanger ({len(hx)})")
        return mask
    if "positions" in spec:
        mask = np.zeros(len(hx), int)
        for xb in spec["positions"]:
            hit = np.flatnonzero(np.abs(hx - float(xb)) < 1e-6)
            if len(hit) == 0:
                raise InputError(f"{section}: restraint at {xb:g} is not at a hanger")
            mask[hit[0]] = 1
        return mask
    return _count_mask(hx, length, parse_number(section, "count", spec["count"], positive=False, integer=True))


def resolve(desc: dict) -> dict:
    """Validated geometry: hanger positions and restraint masks of every line."""
    check_keys("description", desc, {"name", "description", "pipe", "main_line", "branches", "trapezes", "_base_dir"})
    defaults = read_yaml(DEFAULTS_DIR / "system.yaml")
    pipe = deep_merge(defaults["pipe"], desc.get("pipe") or {})
    check_keys("pipe", pipe, set(defaults["pipe"]))
    m = desc["main_line"]
    check_keys("main_line", m, {"length", "n_pipes", "hangers", "braces"})
    L = parse_number("main_line", "length", m.get("length"))
    hx = _hanger_positions("main_line.hangers", m.get("hangers") or {}, L, defaults["hangers"])
    braces = m.get("braces") or {}
    check_keys("main_line.braces", braces, {"transverse", "longitudinal"})
    main = dict(length=L, n_pipes=parse_number("main_line", "n_pipes", m.get("n_pipes"), integer=True),
                hangers=hx,
                trans=_brace_mask("main_line.braces.transverse", braces.get("transverse"), hx, L),
                long=_brace_mask("main_line.braces.longitudinal", braces.get("longitudinal"), hx, L))

    if not desc.get("branches"):
        raise InputError("At least one branch is required")
    branches = []
    for k, b in enumerate(desc["branches"], 1):
        s = f"branches[{k}]"
        check_keys(s, b, {"x", "side", "length", "n_pipes", "hangers", "braces"})
        if str(b.get("side")) not in ("+y", "-y"):
            raise InputError(f"{s}: side must be +y or -y")
        x = parse_number(s, "x", b.get("x"), positive=False)
        if not 0 < x <= L:
            raise InputError(f"{s}: x must lie on the main line (0 < x <= {L:g})")
        Lb = parse_number(s, "length", b.get("length"))
        hy = _hanger_positions(f"{s}.hangers", b.get("hangers") or m.get("hangers") or {}, Lb, defaults["hangers"])
        bb = b.get("braces") or {}
        check_keys(f"{s}.braces", bb, {"transverse", "longitudinal"})
        branches.append(dict(x=x, sign=1.0 if b["side"] == "+y" else -1.0, length=Lb,
                             n_pipes=parse_number(s, "n_pipes", b.get("n_pipes", main["n_pipes"]), integer=True),
                             hangers=hy,
                             trans=_brace_mask(f"{s}.braces.transverse", bb.get("transverse"), hy, Lb),
                             long=_brace_mask(f"{s}.braces.longitudinal", bb.get("longitudinal"), hy, Lb)))

    base = Path(desc.get("_base_dir", "."))
    tz = desc.get("trapezes") or {}
    trap = {k: load_trapeze(base / tz[k] if tz.get(k, "default") != "default" else "default", k)
            for k in ("transverse", "longitudinal")}
    return dict(name=desc["name"], description=desc.get("description", ""), pipe=pipe, main=main,
                branches=branches, trapezes=trap)


# ------------------------------------------------------------------------------------------ model
class _Recorder:
    """Stand-in for openseespy.opensees: records the commands."""

    def __init__(self):
        self.commands = []

    def __getattr__(self, name):
        def call(*args):
            self.commands.append([name, [_plain(a) for a in args]])
        return call


def _plain(v):
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    return v


def build_model3d(desc: dict) -> dict:
    """The 3D model of a system description (load_description) in the format of inputs/models3d."""
    r = resolve(desc)
    op = _Recorder()
    p, main, branches = r["pipe"], r["main"], r["branches"]
    L = main["length"]
    rho_f = p["density"] / 7.8 if p.get("fluid_density") is None else p["fluid_density"]

    op.wipe()
    op.model('basic', '-ndm', 3, '-ndf', 6)

    # ---- hanger (fixed) nodes: main line, then branches
    main_hangers, branch_hangers = [], []
    tag = 1
    for x in main["hangers"]:
        op.node(tag, float(x), 0.0, Z_HANGER)
        op.fix(tag, 1, 1, 1, 1, 1, 1)
        main_hangers.append(tag)
        tag += 1
    for b in branches:
        tags = []
        for d in b["hangers"]:
            op.node(tag, b["x"], b["sign"] * float(d), Z_HANGER)
            op.fix(tag, 1, 1, 1, 1, 1, 1)
            tags.append(tag)
            tag += 1
        branch_hangers.append(tags)

    # ---- masses: (pipe + fluid) per unit length times the mass factor, shared equally by the hangers
    def node_mass(n_pipes, length, n_hangers):
        A, Aw, _, _ = _section(n_pipes, p["outer_diameter"], p["inner_diameter"])
        return p["mass_factor"] * (p["density"] * A + rho_f * Aw) * length / n_hangers

    m_main = node_mass(main["n_pipes"], L, len(main_hangers))
    m_branch = [node_mass(b["n_pipes"], b["length"], len(h)) for b, h in zip(branches, branch_hangers)]

    # ---- main-line pipe nodes: start, hangers, tees, end
    main_nodes, coords, mass_of = [], {}, {}

    def pipe_node(t, x, y, m=0.0):
        if m > 0:
            op.node(t, x, y, 0.0, '-mass', m, m, 0, 0, 0, 0)
        else:
            op.node(t, x, y, 0.0)
        coords[t], mass_of[t] = (x, y), m

    pipe_node(1001, 0.0, 0.0)
    main_nodes.append(1001)
    for h, x in zip(main_hangers, main["hangers"]):
        pipe_node(100 + h, float(x), 0.0, m_main)
        main_nodes.append(100 + h)
    next_tag = 1002
    tees = sorted({b["x"] for b in branches})
    for xt in tees:       # one tee per position; an existing main node at that position is the tee
        if any(abs(coords[t][0] - xt) < 1e-9 for t in main_nodes):
            continue
        pipe_node(next_tag, xt, 0.0)
        main_nodes.append(next_tag)
        next_tag += 1
    if tees[-1] < L:
        pipe_node(next_tag, L, 0.0)
        main_nodes.append(next_tag)
        next_tag += 1
    tee_node = {xt: next(t for t in main_nodes if abs(coords[t][0] - xt) < 1e-9) for xt in tees}

    # ---- branch pipe nodes: hangers, then the free end
    branch_nodes = []
    for j, (b, hs) in enumerate(zip(branches, branch_hangers)):
        tags = []
        for h, d in zip(hs, b["hangers"]):
            pipe_node(100 + h, b["x"], b["sign"] * float(d), m_branch[j])
            tags.append(100 + h)
        pipe_node(next_tag + j, b["x"], b["sign"] * b["length"])
        tags.append(next_tag + j)
        branch_nodes.append(tags)

    # ---- spring nodes: two rows (2xx under the hangers, 3xx above the pipe)
    all_hangers = main_hangers + [h for hs in branch_hangers for h in hs]
    for s in (1, 2):
        for h in all_hangers:
            x, y = coords[100 + h]
            op.node(100 * s + 100 + h, x, y, Z_SPRING)

    # ---- pipe-joint nodes (two coincident sets, 8xx and 9xx, joined by the joint springs)
    main_joints = sorted({k * JOINT_SPACING for k in range(1, int(L // JOINT_SPACING) + 1)} | set(tees))
    branch_joints = [[k * JOINT_SPACING for k in range(1, int(b["length"] // JOINT_SPACING) + 1)
                      if k * JOINT_SPACING < b["length"]] for b in branches]
    joint_sets = []          # per set: main tags, branch tags (list per branch)
    for s in (0, 1):
        t = (s + 1) * 100 + 701
        mt = []
        for xc in main_joints:
            for x in (xc - JOINT_HALF, xc + JOINT_HALF):
                if x <= L:
                    op.node(t, x, 0.0, 0.0)
                    coords[t] = (x, 0.0)
                    mt.append(t)
                    t += 1
        bt = []
        for b, joints in zip(branches, branch_joints):
            tags = []
            for d in [BRANCH_START] + [v for c in joints for v in (c - JOINT_HALF, c + JOINT_HALF)]:
                if d <= b["length"]:
                    op.node(t, b["x"], b["sign"] * d, 0.0)
                    coords[t] = (b["x"], b["sign"] * d)
                    tags.append(t)
                    t += 1
            bt.append(tags)
        joint_sets.append((mt, bt))

    # ---- transformation and materials
    op.geomTransf('Linear', 1, 0, 0, 1)
    op.uniaxialMaterial('Pinching4', 10, *r["trapezes"]["longitudinal"].opensees_args())
    op.uniaxialMaterial('Pinching4', 20, *r["trapezes"]["transverse"].opensees_args())
    op.uniaxialMaterial('MinMax', 1, 10, '-min', -CAP_LONG, '-max', CAP_LONG)
    op.uniaxialMaterial('MinMax', 2, 20, '-min', -CAP_TRANS, '-max', CAP_TRANS)
    op.uniaxialMaterial('Elastic', 3, 0.01)
    op.uniaxialMaterial('Elastic', 4, 10e12)
    op.uniaxialMaterial('ElasticBilin', 33, 141.5, 0.01, 12.8)
    op.uniaxialMaterial('MinMax', 88, 33, '-min', -125, '-max', 125)
    op.uniaxialMaterial('Hysteretic', 7, 7300000, 0.00065, 17600000, 0.0045, 21000000, 0.01,
                        -7300000, -0.00065, -17600000, -0.0045, -21000000, -0.01, 0.8, 0.05, 0, 0, 0.1)
    op.uniaxialMaterial('Parallel', 8, 7, 7, 7)

    def section(n_pipes):
        A, _, J, I = _section(n_pipes, p["outer_diameter"], p["inner_diameter"])
        return [A, p["elastic_modulus"], p["shear_modulus"], J, I, I]

    # ---- pipe elements: main line (sorted by x), then each branch (sorted by distance from the tee)
    ele = [1]

    def chain(nodes, key, sec, joint_tags):
        nodes = sorted(nodes, key=key)
        for a, b in zip(nodes[:-1], nodes[1:]):
            if abs(key(b) - key(a)) < 1e-9:
                if not (a in joint_tags and b in joint_tags):
                    raise InputError(f"Pipe nodes {a} and {b} coincide at {coords[a]}")
                continue
            op.element('elasticBeamColumn', ele[0], a, b, *sec, 1)
            ele[0] += 1

    main_joint_tags = set(joint_sets[0][0]) | set(joint_sets[1][0])
    chain(main_nodes + joint_sets[0][0] + joint_sets[1][0], lambda t: coords[t][0], section(main["n_pipes"]),
          main_joint_tags)
    for j, b in enumerate(branches):
        jt = joint_sets[0][1][j] + joint_sets[1][1][j]
        sec = section(b["n_pipes"])
        op.element('elasticBeamColumn', ele[0], tee_node[b["x"]], joint_sets[0][1][j][0], *sec, 1)
        ele[0] += 1
        chain(branch_nodes[j] + jt, lambda t: abs(coords[t][1]), sec, set(jt))

    # ---- rigid links: fixed hanger node -> spring row 2xx, spring row 3xx -> pipe node
    for h in all_hangers:
        op.rigidLink('beam', h, 200 + h)
    for h in all_hangers:
        op.rigidLink('beam', 300 + h, 100 + h)

    # ---- hanger springs: trapeze (1 longitudinal, 2 transverse) or soft (3) in x and y; rigid otherwise
    braced = {"X": set(), "Y": set()}
    t = 501

    def spring(h, mat_x, mat_y):
        nonlocal t
        op.element('zeroLength', t, 200 + h, 300 + h, '-mat', mat_x, mat_y, 4, 4, 4, 4, '-dir', 1, 2, 3, 4, 5, 6)
        t += 1
        if mat_x != 3:
            braced["X"].add(100 + h)
        if mat_y != 3:
            braced["Y"].add(100 + h)

    for k, h in enumerate(main_hangers):          # main line: longitudinal = x, transverse = y
        spring(h, 1 if main["long"][k] else 3, 2 if main["trans"][k] else 3)
    for b, hs in zip(branches, branch_hangers):   # branches: transverse = x, longitudinal = y
        for k, h in enumerate(hs):
            spring(h, 2 if b["trans"][k] else 3, 1 if b["long"][k] else 3)

    # ---- joint springs (rotational hinge about z) + equal translations and other rotations
    t = 701
    pairs = list(zip(joint_sets[0][0], joint_sets[1][0]))
    for j in range(len(branches)):
        pairs += list(zip(joint_sets[0][1][j], joint_sets[1][1][j]))
    for a, b in pairs:
        op.element('zeroLength', t, a, b, '-mat', 4, 4, 4, 4, 4, 8, '-dir', 1, 2, 3, 4, 5, 6)
        op.equalDOF(a, b, 1, 2, 3, 4, 5)
        t += 1

    # ---- gravity: weight of each node's own mass
    op.timeSeries('Constant', 1)
    op.pattern('Plain', 1, 1)
    for n in main_nodes + [n for ns in branch_nodes for n in ns]:
        if mass_of[n] > 0:
            op.load(n, 0.0, 0.0, -mass_of[n] * G_ACC, 0.0, 0.0, 0.0)

    # ---- recorded nodes and restrained nodes
    nodes = sorted(main_nodes, key=lambda n: coords[n][0])
    for ns in branch_nodes:
        nodes += sorted(ns, key=lambda n: coords[n][1])
    groups = {"X": nodes, "Y": list(nodes)}
    for d in "XY":
        groups[d + "t"] = [n for n in nodes if n in braced[d]]
        groups[d + "p"] = [n for n in nodes if n not in braced[d]]

    return {"name": r["name"], "description": r["description"] or f"3D model of {r['name']} (automatized workflow)",
            "units": "N, mm, s", "rom_systems": {}, "node_groups": groups, "commands": op.commands}


def write_model3d(model: dict, path: str | Path) -> None:
    Path(path).write_text(json.dumps(model, indent=None, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Build the 3D model of a system description")
    ap.add_argument("description", help="system description (YAML)")
    ap.add_argument("-o", "--output", help="output JSON (default: <name>.json next to the description)")
    a = ap.parse_args()
    m = build_model3d(load_description(a.description))
    out = Path(a.output) if a.output else Path(a.description).with_name(f"{m['name']}.json")
    write_model3d(m, out)
    g = m["node_groups"]
    print(f"{out}: {sum(c == 'node' for c, _ in m['commands'])} nodes, "
          f"{sum(c == 'element' for c, _ in m['commands'])} elements; recorded nodes {len(g['X'])}, "
          f"restrained in x {len(g['Xt'])}, in y {len(g['Yt'])}")
