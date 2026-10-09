"""Verification with the full 3D models of the paper (``inputs/models3d/<name>.json``).

The model files hold the OpenSees commands of the paper's 3D models, unchanged (converted by
``validation/convert_3d_models.py``), and the test types of each script's fallback sequence. The analysis
follows ``3D_models/*_biron.py``: gravity, mass-proportional damping at the first (elastic) mode,
bidirectional floor motions (x and y components applied together), Newmark average acceleration in chunks
that stop when a braced node exceeds the collapse displacement, and the same fallback sequence.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import openseespy.opensees as op

from . import DEFAULTS_DIR, INPUTS_DIR
from .inputs import InputError, check_keys, deep_merge, parse_number, parse_rayleigh, read_yaml
from .motions import FloorMotion

MODELS3D_DIR = INPUTS_DIR / "models3d"


@dataclass
class Model3D:
    name: str
    description: str
    commands: list
    node_groups: dict[str, list[int]]
    rom_systems: dict[str, str]
    fallback_tests: list[list[str]] | None = None   # per fallback attempt: [test during it, main test after it]

    def node_coords(self) -> dict[int, tuple[float, float, float]]:
        return {int(a[0]): tuple(float(v) for v in a[1:4]) for c, a in self.commands if c == "node"}


def list_models3d() -> list[str]:
    return sorted(p.stem for p in MODELS3D_DIR.glob("*.json"))


def load_model3d(name: str) -> Model3D:
    path = MODELS3D_DIR / f"{name}.json"
    if not path.exists():
        raise InputError(f"No 3D model '{name}' in inputs/models3d/ (available: {', '.join(list_models3d())})")
    d = json.loads(path.read_text())
    return Model3D(d["name"], d.get("description", ""), d["commands"], d["node_groups"], d.get("rom_systems", {}),
                   d.get("fallback_tests"))


def model3d_for_system(system_name: str) -> str | None:
    for name in list_models3d():
        if system_name in load_model3d(name).rom_systems.values():
            return name
    return None


@dataclass
class Verification3DSettings:
    gravity_steps: int
    damping_ratio: float
    rayleigh_modes: tuple[int, ...]
    mass_proportional: float
    current_stiffness: float
    committed_stiffness: float
    initial_stiffness: float
    system: str
    test: str
    tolerance: float
    max_iterations: int
    fallback_test: str | None           # None: the test types of the model's own script
    fallback_tolerance: float
    fallback_max_iterations: int
    collapse_displacement: float | None  # largest braced-node displacement (mm) that stops the analysis
    check_steps: int
    record_dt: float | None              # time step of the displacement recorders (None: every step)

    @classmethod
    def from_dict(cls, data: dict | None = None) -> "Verification3DSettings":
        d = deep_merge(read_yaml(DEFAULTS_DIR / "settings.yaml").get("verification_3d", {}), data or {})
        s = "verification_3d"
        check_keys(s, d, {"gravity_steps", "damping_ratio", "rayleigh_modes", "rayleigh", "system", "test",
                          "tolerance", "max_iterations", "fallback_test", "fallback_tolerance",
                          "fallback_max_iterations", "collapse_displacement", "check_steps", "record_dt"})
        mass, current, committed, initial = parse_rayleigh(s, d["rayleigh"])
        modes = tuple(int(v) for v in d["rayleigh_modes"])
        if len(modes) not in (1, 2) or min(modes) < 1:
            raise InputError("'verification_3d.rayleigh_modes' must be one or two mode numbers, e.g. [1] or [1, 2]")
        return cls(
            gravity_steps=parse_number(s, "gravity_steps", d["gravity_steps"], integer=True),
            damping_ratio=parse_number(s, "damping_ratio", d["damping_ratio"], positive=False),
            rayleigh_modes=modes,
            mass_proportional=mass, current_stiffness=current, committed_stiffness=committed,
            initial_stiffness=initial,
            system=str(d["system"]),
            test=str(d["test"]), tolerance=parse_number(s, "tolerance", d["tolerance"]),
            max_iterations=parse_number(s, "max_iterations", d["max_iterations"], integer=True),
            fallback_test=None if d["fallback_test"] in (None, "script") else str(d["fallback_test"]),
            fallback_tolerance=parse_number(s, "fallback_tolerance", d["fallback_tolerance"]),
            fallback_max_iterations=parse_number(s, "fallback_max_iterations", d["fallback_max_iterations"], integer=True),
            collapse_displacement=None if d["collapse_displacement"] is None
            else parse_number(s, "collapse_displacement", d["collapse_displacement"]),
            check_steps=parse_number(s, "check_steps", d["check_steps"], integer=True),
            record_dt=None if d["record_dt"] is None else parse_number(s, "record_dt", d["record_dt"]),
        )


@dataclass
class Response3D:
    model: str
    record_x: str
    record_y: str
    level: int | None
    periods: list[float]
    time: np.ndarray
    nodes_x: list[int]
    nodes_y: list[int]
    ux: np.ndarray            # (n_time, n_nodes_x) displacement in x relative to the floor (mm)
    uy: np.ndarray            # (n_time, n_nodes_y)
    completed: bool
    end_time: float
    duration: float
    collapsed: bool = False   # stopped because a braced node exceeded the collapse displacement

    @property
    def peak_x(self) -> np.ndarray:
        return np.max(np.abs(self.ux), axis=0) if self.ux.size else np.full(len(self.nodes_x), np.nan)

    @property
    def peak_y(self) -> np.ndarray:
        return np.max(np.abs(self.uy), axis=0) if self.uy.size else np.full(len(self.nodes_y), np.nan)


def build_model3d(model: Model3D) -> None:
    op.wipe()
    op.logFile(os.devnull, "-noEcho")
    for cmd, args in model.commands:
        if cmd == "wipe":
            continue
        getattr(op, cmd)(*args)


def run_3d(model: Model3D, motion_x: FloorMotion, motion_y: FloorMotion, s: Verification3DSettings,
           progress: Callable[[float], None] | None = None) -> Response3D:
    if abs(motion_x.dt - motion_y.dt) > 1e-12:
        raise InputError("The x and y floor motions must have the same time step")
    build_model3d(model)

    # gravity (static, load control)
    op.constraints("Transformation")
    op.numberer("RCM")
    op.system("BandGeneral")
    op.test("NormDispIncr", 1.0e-8, 200)
    op.algorithm("Newton")
    op.integrator("LoadControl", 1.0 / s.gravity_steps)
    op.analysis("Static")
    if op.analyze(s.gravity_steps) != 0:
        raise RuntimeError("Gravity analysis of the 3D model failed")
    op.loadConst("-time", 0.0)

    periods = [2 * math.pi / math.sqrt(lam) for lam in op.eigen(4)]
    xi = s.damping_ratio
    if len(s.rayleigh_modes) == 1:
        # one mode: each term alone gives xi at that mode (2 xi w on the mass, 2 xi / w on the stiffness)
        w = op.eigen(s.rayleigh_modes[0])[s.rayleigh_modes[0] - 1] ** 0.5   # as the scripts
        a_m, a_k = 2 * xi * w, 2 * xi / w
    else:
        i, j = s.rayleigh_modes
        lam = op.eigen(max(i, j))
        wi, wj = math.sqrt(lam[i - 1]), math.sqrt(lam[j - 1])
        a_m, a_k = xi * (2 * wi * wj) / (wi + wj), 2 * xi / (wi + wj)
    op.rayleigh(s.mass_proportional * a_m, s.current_stiffness * a_k, s.initial_stiffness * a_k,
                s.committed_stiffness * a_k)

    dt = motion_x.dt
    npts = min(len(motion_x.acc), len(motion_y.acc))
    t_max = dt * npts
    op.timeSeries("Path", 1000, "-dt", dt, "-values", *motion_x.acc, "-factor", 1, "-prependZero")
    op.timeSeries("Path", 2000, "-dt", dt, "-values", *motion_y.acc, "-factor", 1, "-prependZero")
    op.pattern("UniformExcitation", 100, 1, "-accel", 1000)
    op.pattern("UniformExcitation", 200, 2, "-accel", 2000)

    nodes_x, nodes_y = model.node_groups["X"], model.node_groups["Y"]
    braced_x, braced_y = model.node_groups.get("Xt", []), model.node_groups.get("Yt", [])

    def over_cap() -> bool:
        if s.collapse_displacement is None:
            return False
        return max([abs(op.nodeDisp(n, 1)) for n in braced_x] + [abs(op.nodeDisp(n, 2)) for n in braced_y]
                   + [0.0]) > s.collapse_displacement

    # test types of the four fallback attempts: [during, main test restored after]
    if s.fallback_test is not None or not model.fallback_tests:
        fb = [[s.fallback_test or s.test, s.test]] * 4
    else:
        fb = model.fallback_tests

    with tempfile.TemporaryDirectory() as tmp:
        fx, fy = Path(tmp) / "x.out", Path(tmp) / "y.out"
        rec_dt = [] if s.record_dt is None else ["-dt", s.record_dt]
        op.recorder("Node", "-file", str(fx), "-time", *rec_dt, "-node", *nodes_x, "-dof", 1, "disp")
        op.recorder("Node", "-file", str(fy), "-time", *rec_dt, "-node", *nodes_y, "-dof", 2, "disp")

        op.wipeAnalysis()
        op.integrator("Newmark", 0.5, 0.25)
        op.numberer("RCM")
        op.system(s.system)
        op.constraints("Transformation")
        op.test(s.test, s.tolerance, s.max_iterations)
        op.algorithm("Newton")
        op.analysis("Transient")

        # in chunks of check_steps steps, stopped if a braced node exceeds the collapse displacement
        collapsed = False
        done, ok = 0, 0
        while done < npts and ok == 0:
            n = min(s.check_steps, npts - done)
            ok = op.analyze(n, dt)
            done += n
            if progress:
                progress(min(op.getTime() / t_max, 1.0))
            if ok == 0 and over_cap():
                collapsed = True
                break

        if ok != 0 and not collapsed:
            def attempt(k: int, step: float) -> int:
                op.test(fb[k][0], s.fallback_tolerance, s.fallback_max_iterations)
                res = op.analyze(1, step)
                op.test(fb[k][1], s.tolerance, s.max_iterations)
                op.algorithm("Newton")
                return res

            ok = 0
            control_time = op.getTime()
            last_report = control_time
            while control_time < t_max and ok == 0:
                control_time = op.getTime()
                if over_cap():
                    collapsed = True
                    break
                ok = op.analyze(1, dt)
                if ok != 0:
                    op.algorithm("Newton", "-initial")
                    ok = attempt(0, dt / 2)
                if ok != 0:
                    op.algorithm("Broyden", 50)
                    ok = attempt(1, dt / 2)
                if ok != 0:
                    op.algorithm("NewtonLineSearch")
                    ok = attempt(2, dt / 10)
                if ok != 0:
                    op.algorithm("KrylovNewton")
                    ok = attempt(3, dt / 20)
                if progress and control_time - last_report > 0.02 * t_max:
                    last_report = control_time
                    progress(min(control_time / t_max, 1.0))
        collapsed = collapsed or over_cap()     # also after the last step
        end_time = op.getTime()
        op.wipe()
        ax = np.loadtxt(fx, ndmin=2)
        ay = np.loadtxt(fy, ndmin=2)

    n = min(len(ax), len(ay))
    return Response3D(
        model=model.name, record_x=motion_x.record, record_y=motion_y.record, level=motion_x.level,
        periods=periods, time=ax[:n, 0], nodes_x=list(nodes_x), nodes_y=list(nodes_y),
        ux=ax[:n, 1:], uy=ay[:n, 1:], completed=end_time >= t_max - 0.5 * dt, end_time=end_time, duration=t_max,
        collapsed=collapsed,
    )


def braced_peaks(model: Model3D, response: Response3D) -> dict[str, float]:
    """Largest peak displacement over the braced nodes, per direction: the measure stored in the paper's
    ``Results/*_DispX.txt`` / ``*_DispY.txt`` (where a collapsed run is stored as 100 mm and a run that did
    not reach the end of the motion as 0)."""
    out = {}
    for d, nodes, peaks in (("x", response.nodes_x, response.peak_x), ("y", response.nodes_y, response.peak_y)):
        idx = [nodes.index(n) for n in model.node_groups.get(f"{d.upper()}t", []) if n in nodes]
        out[d] = float(peaks[idx].max()) if idx else float("nan")
    return out
