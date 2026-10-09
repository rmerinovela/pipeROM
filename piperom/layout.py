"""Plan layout of a piping system and its reduction to the equivalent system of each direction.

A layout is one straight main line along x (from x = 0) and straight branches along y framing into it, on
either side (+y or -y), with gravity hangers and trapeze restraints in both directions on every line. The
format is the one of ``automatized_workflow/model3d.py`` (which builds the 3D model from it):

    name: MySystem
    main_line:
      length: 36000
      n_pipes: 3
      hangers: {first: 1000, spacing: 3000}      # or {positions: [...]}; optional end_clearance
      braces:
        transverse: {count: 4}                   # along y; or {mask: [...]} / {positions: [...]}
        longitudinal: {count: 3}                 # along x
    branches:
      - x: 36000                                 # tee position on the main line (0 < x <= length)
        side: +y                                 # +y or -y
        length: 18000
        n_pipes: 3                               # optional, default: main line's
        hangers: {first: 1000, spacing: 3000}    # optional, default: main line's; distances from the tee
        braces:
          transverse: {count: 2}                 # along x
          longitudinal: {count: 2}               # along y
    pipe: {...}                                  # optional, as in the system files
    trapezes: {...}                              # optional, as in the system files

Rules of the reduced-order model (as the paper's 2D models consistent with their 3D models, e.g. M29-M31):

* all branches have the same length, and the branches on each side of the main line are identical (pipes,
  hangers and restraints); the tee positions are free;
* **loading along y**: the analysed line is the main line with its transverse (y) restraints; each tee is
  one branch DOF with the branch's longitudinal (y) restraints; branches sharing a tee (one on each side)
  are one DOF with their pipes and restraints added;
* **loading along x**: the analysed line is a branch line through the main line, at its tee: with branches
  on both sides, the -y branch to the left of the tee and the +y branch to the right (line of twice the
  branch length, tee at the centre); with one side only, the branch to the left of the tee (tee at the end
  of the line). The branches of each side are lumped as identical mains (``n_mains_left`` / ``n_mains_right``,
  split at the tee), and the main line is the single branch DOF at the tee, with its longitudinal (x)
  restraints.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import DEFAULTS_DIR
from .inputs import InputError, PipingSystem, check_keys, evenly_spaced_brace_mask, parse_number, read_yaml

SIDES = ("+y", "-y")


@dataclass
class Line:
    """One line of the layout: hangers (distance from its start: x = 0 for the main line, the tee for a
    branch) and restraint masks per direction of the line."""

    length: float
    n_pipes: int
    hangers: np.ndarray
    transverse: np.ndarray       # restraints across the line
    longitudinal: np.ndarray     # restraints along the line

    def same_as(self, other: "Line") -> bool:
        return (self.length == other.length and self.n_pipes == other.n_pipes
                and np.array_equal(self.hangers, other.hangers)
                and np.array_equal(self.transverse, other.transverse)
                and np.array_equal(self.longitudinal, other.longitudinal))


@dataclass
class Branch:
    x: float
    side: str
    line: Line


@dataclass
class Layout:
    name: str
    description: str
    main: Line
    branches: list[Branch]
    pipe: dict = field(default_factory=dict)        # overrides of inputs/defaults/system.yaml
    trapezes: dict = field(default_factory=dict)
    raw: dict = field(default_factory=dict)         # the description as given (for model3d and saving)

    # ------------------------------------------------------------------ geometry
    def sides(self) -> dict[str, list[Branch]]:
        return {s: [b for b in self.branches if b.side == s] for s in SIDES}

    def branch_line(self, side: str) -> Line | None:
        bs = self.sides()[side]
        return bs[0].line if bs else None

    # ------------------------------------------------------------------ reduction
    def system_y(self) -> PipingSystem:
        """Equivalent system for loading along y (the main line analysed)."""
        sides = self.sides()
        tees = sorted({b.x for b in sides["+y"]}) + sorted({b.x for b in sides["-y"]} - {b.x for b in sides["+y"]})
        branches = []
        for x in tees:
            at = [b for b in self.branches if b.x == x]
            n_braces = int(sum(b.line.longitudinal.sum() for b in at))
            if n_braces == 0:
                raise InputError(f"Loading along y: the branch at x = {x:g} mm has no restraint along y "
                                 "(longitudinal restraint of the branch); the reduced-order model needs at least one")
            branches.append({"x": _num(x), "length": _num(at[0].line.length),
                             "n_pipes": sum(b.line.n_pipes for b in at), "n_braces": n_braces})
        return self._system(f"{self.name}y", "y", self.main.length, self.main.n_pipes, self.main.hangers,
                            self.main.transverse, branches, {})

    def system_x(self) -> PipingSystem:
        """Equivalent system for loading along x (a branch line analysed, branches lumped as mains)."""
        sides = self.sides()
        neg, pos = self.branch_line("-y"), self.branch_line("+y")
        Lb = (neg or pos).length
        n_pipes = (neg or pos).n_pipes
        hangers, mask = [], []
        if neg is not None and pos is not None:
            length, tee = 2 * Lb, Lb
            for d, t in sorted(zip(neg.hangers, neg.transverse), reverse=True):
                hangers.append(tee - d); mask.append(int(t))
            for d, t in zip(pos.hangers, pos.transverse):
                hangers.append(tee + d); mask.append(int(t))
            mains = {"n_mains_left": len(sides["-y"]), "n_mains_right": len(sides["+y"]), "x_center": _num(tee)}
        else:
            line = neg or pos
            length, tee = Lb, Lb
            for d, t in sorted(zip(line.hangers, line.transverse), reverse=True):
                hangers.append(tee - d); mask.append(int(t))
            mains = {"n_mains": len(self.branches)}
        if self.main.longitudinal.sum() == 0:
            raise InputError("Loading along x: the main line has no restraint along x (longitudinal restraint of the "
                             "main line); the reduced-order model needs at least one")
        branches = [{"x": _num(tee), "length": _num(self.main.length), "n_pipes": self.main.n_pipes,
                     "n_braces": int(self.main.longitudinal.sum())}]
        return self._system(f"{self.name}x", "x", length, n_pipes, np.array(hangers), np.array(mask), branches, mains)

    def systems(self) -> dict[str, PipingSystem]:
        return {"x": self.system_x(), "y": self.system_y()}

    def _system(self, name, direction, length, n_pipes, hangers, mask, branches, mains) -> PipingSystem:
        data = {"name": name,
                "description": f"Loading along {direction} of layout {self.name} (piperom.layout)",
                "main_line": {"length": _num(length), "n_pipes": int(n_pipes), **mains},
                "hangers": {"positions": [_num(h) for h in hangers]},
                "braces": {"mask": [int(v) for v in mask]},
                "branches": branches}
        if self.pipe:
            data["pipe"] = dict(self.pipe)
        if self.trapezes:
            data["trapezes"] = dict(self.trapezes)
        return PipingSystem.from_dict(data, base_dir=Path(self.raw.get("_base_dir", ".")))


# ----------------------------------------------------------------------------- parsing
def _num(v):
    v = float(v)
    return int(v) if v.is_integer() else v


def _hanger_positions(section: str, spec: dict, length: float, end_clearance: float) -> np.ndarray:
    """Hanger positions along a line (distance from its start), with the rules of the system files."""
    check_keys(section, spec, {"first", "spacing", "end_clearance", "positions"})
    if spec.get("positions") is not None:
        x = np.round(np.array(spec["positions"], float), 6)
        if len(x) == 0 or np.any(np.diff(x) <= 0) or x[0] <= 0 or x[-1] >= length:
            raise InputError(f"{section}: positions must be increasing and strictly inside the line")
        return x
    x0 = parse_number(section, "first", spec.get("first"))
    s = parse_number(section, "spacing", spec.get("spacing"))
    c = parse_number(section, "end_clearance", spec.get("end_clearance", end_clearance), positive=False)
    if x0 > length - c:
        raise InputError(f"{section}: the first hanger must be at least {c:g} mm from the end")
    n_max = int((length - c - x0) // s)
    return np.round(np.array([x0 + i * s for i in range(n_max + 1) if x0 + i * s <= length - c + 1e-6]), 6)


def _brace_mask(section: str, spec: dict | None, hx: np.ndarray, length: float) -> np.ndarray:
    spec = spec or {"count": 0}
    check_keys(section, spec, {"mask", "positions", "count"})
    if sum(k in spec for k in ("mask", "positions", "count")) != 1:
        raise InputError(f"{section}: give exactly one of mask, positions, count")
    if "mask" in spec:
        mask = np.array(spec["mask"], int)
        if len(mask) != len(hx) or not np.all((mask == 0) | (mask == 1)):
            raise InputError(f"{section}: the mask needs one 0/1 per hanger ({len(hx)})")
        return mask
    if "positions" in spec:
        mask = np.zeros(len(hx), int)
        for xb in spec["positions"]:
            hit = np.flatnonzero(np.abs(hx - float(xb)) < 1e-6)
            if len(hit) == 0:
                raise InputError(f"{section}: restraint at {float(xb):g} is not at a hanger")
            mask[hit[0]] = 1
        return mask
    return evenly_spaced_brace_mask(hx, length, parse_number(section, "count", spec["count"], positive=False,
                                                             integer=True))


def _line(section: str, d: dict, length: float, n_pipes: int, hanger_spec: dict, end_clearance: float) -> Line:
    hx = _hanger_positions(f"{section}.hangers", hanger_spec, length, end_clearance)
    braces = d.get("braces") or {}
    check_keys(f"{section}.braces", braces, {"transverse", "longitudinal"})
    return Line(length=float(length), n_pipes=n_pipes, hangers=hx,
                transverse=_brace_mask(f"{section}.braces.transverse", braces.get("transverse"), hx, length),
                longitudinal=_brace_mask(f"{section}.braces.longitudinal", braces.get("longitudinal"), hx, length))


def layout_from_dict(desc: dict, base_dir: Path | None = None) -> Layout:
    check_keys("layout", desc, {"name", "description", "pipe", "main_line", "branches", "trapezes", "_base_dir"})
    defaults = read_yaml(DEFAULTS_DIR / "system.yaml")
    clearance = float(defaults["hangers"]["end_clearance"])
    m = desc.get("main_line") or {}
    check_keys("main_line", m, {"length", "n_pipes", "hangers", "braces"})
    L = parse_number("main_line", "length", m.get("length"))
    n_main = parse_number("main_line", "n_pipes", m.get("n_pipes"), integer=True)
    main = _line("main_line", m, L, n_main, m.get("hangers") or {}, clearance)

    if not desc.get("branches"):
        raise InputError("At least one branch is required")
    branches = []
    for k, b in enumerate(desc["branches"], 1):
        s = f"branches[{k}]"
        check_keys(s, b, {"x", "side", "length", "n_pipes", "hangers", "braces"})
        side = str(b.get("side"))
        if side not in SIDES:
            raise InputError(f"{s}: side must be +y or -y")
        x = parse_number(s, "x", b.get("x"), positive=False)
        if not 0 < x <= L:
            raise InputError(f"{s}: the tee must lie on the main line (0 < x <= {L:g})")
        Lb = parse_number(s, "length", b.get("length"))
        n_b = parse_number(s, "n_pipes", b.get("n_pipes", n_main), integer=True)
        branches.append(Branch(round(float(x), 6), side, _line(s, b, Lb, n_b, b.get("hangers") or m.get("hangers") or {},
                                                               clearance)))

    lengths = {b.line.length for b in branches}
    if len(lengths) > 1:
        raise InputError(f"All branches must have the same length (found {', '.join(f'{v:g}' for v in sorted(lengths))} mm)")
    for side in SIDES:
        bs = [b for b in branches if b.side == side]
        for b in bs[1:]:
            if not b.line.same_as(bs[0].line):
                raise InputError(f"The branches on the {side} side must be identical (pipes, hangers and restraints): "
                                 f"the branch at x = {b.x:g} differs from the one at x = {bs[0].x:g}")
    if {b.line.n_pipes for b in branches} != {branches[0].line.n_pipes}:
        raise InputError("Branches on the two sides must have the same number of pipes")
    seen = set()
    for b in branches:
        if (b.x, b.side) in seen:
            raise InputError(f"Two branches on the {b.side} side at the same tee (x = {b.x:g})")
        seen.add((b.x, b.side))

    return Layout(name=str(desc.get("name") or "layout"), description=str(desc.get("description") or ""),
                  main=main, branches=branches, pipe=dict(desc.get("pipe") or {}),
                  trapezes=dict(desc.get("trapezes") or {}),
                  raw={**desc, "_base_dir": str(base_dir or desc.get("_base_dir") or ".")})


def load_layout(path: str | Path) -> Layout:
    path = Path(path)
    return layout_from_dict(read_yaml(path), base_dir=path.parent)
