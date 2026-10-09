"""Nonlinear time-history analysis of the equivalent SDOF under a floor acceleration history.

Port of ``code implementation for paper/Pushover_SDOF/<tag>_SDOF_NLTHA.py``: one zeroLength element with
the parallel Pinching4 springs of the supports, each wrapped in a MinMax at its restraint cap; mass =
effective mass; mass-proportional damping at the initial (elastic) frequency; Newmark average
acceleration, run in chunks so that the analysis stops when the largest braced-node displacement
Gamma * max(phi) * u exceeds the collapse displacement; and the scripts' fallback sequence.
"""

from __future__ import annotations

import math
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import openseespy.opensees as op

from .inputs import check_keys, deep_merge, parse_number, parse_rayleigh, read_yaml
from . import DEFAULTS_DIR
from .motions import FloorMotion
from .trapeze import Pinching4

RIGID = 10e12


@dataclass
class SDOFModel:
    """What the time-history analysis needs: mass and the springs acting in parallel."""

    name: str
    mass: float
    springs: list[Pinching4]
    gamma: float | None = None
    support_phi: list[float] = field(default_factory=list)      # same order as springs
    support_labels: list[str] = field(default_factory=list)
    support_kinds: list[str] = field(default_factory=list)      # "longitudinal" / "transverse"

    @classmethod
    def from_parameters(cls, p) -> "SDOFModel":
        """Longitudinal springs first, then transverse (order of the SDOF scripts)."""
        sup = [s for s in p.supports if s.kind == "longitudinal"] + [s for s in p.supports if s.kind == "transverse"]
        labels = [f"{s.kind} x={s.x:g}" + ("" if s.branch is None else f" (branch {s.branch + 1})") for s in sup]
        return cls(name=p.name, mass=p.effective_mass, springs=[s.spring for s in sup], gamma=p.gamma,
                   support_phi=[s.phi for s in sup], support_labels=labels, support_kinds=[s.kind for s in sup])

    @property
    def braced_factor(self) -> float | None:
        """Largest braced-node displacement per unit SDOF displacement, Gamma * max(phi)."""
        if self.gamma is None or not self.support_phi:
            return None
        return self.gamma * max(self.support_phi)


@dataclass
class TimeHistorySettings:
    damping_ratio: float
    mass_proportional: float
    current_stiffness: float
    committed_stiffness: float
    initial_stiffness: float
    test: str
    tolerance: float
    max_iterations: int
    fallback_tolerance: float
    fallback_max_iterations: int
    cap_longitudinal: float | None        # restraint caps (mm at the support); None = no cap
    cap_transverse: float | None
    collapse_displacement: float | None   # largest braced-node displacement (mm) that stops the analysis
    check_steps: int

    @classmethod
    def from_dict(cls, data: dict | None = None) -> "TimeHistorySettings":
        d = deep_merge(read_yaml(DEFAULTS_DIR / "settings.yaml").get("sdof_time_history", {}), data or {})
        s = "sdof_time_history"
        check_keys(s, d, {"damping_ratio", "rayleigh", "test", "tolerance", "max_iterations", "fallback_tolerance",
                          "fallback_max_iterations", "restraint_caps", "collapse_displacement", "check_steps"})
        mass, current, committed, initial = parse_rayleigh(s, d["rayleigh"])
        caps = d["restraint_caps"] or {}
        check_keys(f"{s}.restraint_caps", caps, {"longitudinal", "transverse"})

        def optional(section, key, value):
            return None if value is None else parse_number(section, key, value)

        return cls(
            damping_ratio=parse_number(s, "damping_ratio", d["damping_ratio"], positive=False),
            mass_proportional=mass, current_stiffness=current, committed_stiffness=committed,
            initial_stiffness=initial,
            test=str(d["test"]), tolerance=parse_number(s, "tolerance", d["tolerance"]),
            max_iterations=parse_number(s, "max_iterations", d["max_iterations"], integer=True),
            fallback_tolerance=parse_number(s, "fallback_tolerance", d["fallback_tolerance"]),
            fallback_max_iterations=parse_number(s, "fallback_max_iterations", d["fallback_max_iterations"], integer=True),
            cap_longitudinal=optional(f"{s}.restraint_caps", "longitudinal", caps.get("longitudinal")),
            cap_transverse=optional(f"{s}.restraint_caps", "transverse", caps.get("transverse")),
            collapse_displacement=optional(s, "collapse_displacement", d["collapse_displacement"]),
            check_steps=parse_number(s, "check_steps", d["check_steps"], integer=True),
        )

    def to_dict(self) -> dict:
        return {"damping_ratio": self.damping_ratio,
                "rayleigh": {"mass": self.mass_proportional, "current_stiffness": self.current_stiffness,
                             "committed_stiffness": self.committed_stiffness,
                             "initial_stiffness": self.initial_stiffness},
                "test": self.test, "tolerance": self.tolerance, "max_iterations": self.max_iterations,
                "fallback_tolerance": self.fallback_tolerance,
                "fallback_max_iterations": self.fallback_max_iterations,
                "restraint_caps": {"longitudinal": self.cap_longitudinal, "transverse": self.cap_transverse},
                "collapse_displacement": self.collapse_displacement, "check_steps": self.check_steps}


@dataclass
class SDOFResponse:
    motion: str
    record: str
    level: int | None
    period: float
    time: np.ndarray
    u: np.ndarray               # SDOF displacement relative to the floor (mm)
    force: np.ndarray           # spring force (N)
    completed: bool
    end_time: float
    duration: float
    collapsed: bool = False     # stopped because a braced node exceeded the collapse displacement
    peak_u: float = float("nan")   # peak |u| (kept when the histories are dropped)

    def __post_init__(self):
        if len(self.u):
            self.peak_u = float(np.max(np.abs(self.u)))


def _cap(model: SDOFModel, k: int, s: TimeHistorySettings | None) -> float | None:
    """Cap of spring k in SDOF displacement: cap at the support / (Gamma * phi)."""
    if s is None or model.gamma is None or k >= len(model.support_kinds):
        return None
    cap = {"longitudinal": s.cap_longitudinal, "transverse": s.cap_transverse}[model.support_kinds[k]]
    return None if cap is None else cap / (model.gamma * model.support_phi[k])


def _build(model: SDOFModel, s: TimeHistorySettings | None = None) -> None:
    op.wipe()
    op.logFile(os.devnull, "-noEcho")
    op.model("basic", "-ndm", 2, "-ndf", 3)
    op.node(1, 0.0, 0.0)
    op.node(2, 0.0, 0.0, "-mass", model.mass, model.mass, 0)
    op.fix(1, 1, 1, 1)
    tags = []
    for k, spring in enumerate(model.springs):
        tag = 100 + k
        op.uniaxialMaterial("Pinching4", tag, *spring.opensees_args())
        cap = _cap(model, k, s)
        if cap is not None:
            op.uniaxialMaterial("MinMax", 2000 + tag, tag, "-min", -cap, "-max", cap)
            tag += 2000
        tags.append(tag)
    op.uniaxialMaterial("Parallel", 1000, *tags)
    op.uniaxialMaterial("Elastic", 4, RIGID)
    op.element("zeroLength", 1, 1, 2, "-mat", 1000, 4, 4, "-dir", 1, 2, 3)


def sdof_period(model: SDOFModel) -> float:
    _build(model)
    return 2 * math.pi / math.sqrt(op.eigen(1)[0])


def run_sdof_time_history(model: SDOFModel, motion: FloorMotion, s: TimeHistorySettings) -> SDOFResponse:
    _build(model, s)
    omega = op.eigen(1)[0] ** 0.5       # as the scripts (not math.sqrt), for identical damping
    xi = s.damping_ratio
    op.rayleigh(s.mass_proportional * 2 * xi * omega, s.current_stiffness * 2 * xi / omega,
                s.initial_stiffness * 2 * xi / omega, s.committed_stiffness * 2 * xi / omega)

    dt, npts = motion.dt, len(motion.acc)
    t_max = dt * npts
    op.timeSeries("Path", 1000, "-dt", dt, "-values", *motion.acc, "-factor", 1, "-prependZero")
    op.pattern("UniformExcitation", 100, 1, "-accel", 1000)
    factor = model.braced_factor

    def over_cap() -> bool:
        if s.collapse_displacement is None or factor is None:
            return False
        return factor * abs(op.nodeDisp(2, 1)) > s.collapse_displacement

    def main_test():
        op.test(s.test, s.tolerance, s.max_iterations)

    def fallback_test():
        op.test(s.test, s.fallback_tolerance, s.fallback_max_iterations)

    with tempfile.TemporaryDirectory() as tmp:
        f_u, f_f = Path(tmp) / "u.out", Path(tmp) / "f.out"
        op.recorder("Node", "-file", str(f_u), "-time", "-node", 2, "-dof", 1, "disp")
        op.recorder("Element", "-file", str(f_f), "-ele", 1, "force")

        op.wipeAnalysis()
        op.integrator("Newmark", 0.5, 0.25)
        op.numberer("RCM")
        op.system("FullGeneral")
        op.constraints("Transformation")
        main_test()
        op.algorithm("Newton")
        op.analysis("Transient")

        # in chunks of check_steps steps, stopped if a braced node exceeds the collapse displacement
        collapsed = False
        done, ok = 0, 0
        while done < npts and ok == 0:
            n = min(s.check_steps, npts - done)
            ok = op.analyze(n, dt)
            done += n
            if ok == 0 and over_cap():
                collapsed = True
                break

        if ok != 0 and not collapsed:
            ok = 0
            control_time = op.getTime()
            while control_time < t_max and ok == 0:
                control_time = op.getTime()
                if over_cap():
                    collapsed = True
                    break
                ok = op.analyze(1, dt)
                if ok != 0:
                    fallback_test()
                    op.algorithm("Newton", "-initial")
                    ok = op.analyze(1, dt / 2)
                    main_test()
                    op.algorithm("Newton")
                if ok != 0:
                    op.algorithm("Broyden", 50)
                    fallback_test()
                    ok = op.analyze(1, dt / 2)
                    main_test()
                    op.algorithm("Newton")
                if ok != 0:
                    op.algorithm("NewtonLineSearch")
                    fallback_test()
                    ok = op.analyze(1, dt / 10)
                    main_test()
                    op.algorithm("Newton")
                if ok != 0:
                    op.algorithm("KrylovNewton")
                    fallback_test()
                    ok = op.analyze(1, dt / 20)
                    main_test()
                    op.algorithm("Newton")
        collapsed = collapsed or over_cap()     # also after the last step
        end_time = op.getTime()
        op.wipe()
        u = np.loadtxt(f_u, ndmin=2)
        frc = np.loadtxt(f_f, ndmin=2)

    n = min(len(u), len(frc))
    return SDOFResponse(
        motion=motion.set_name, record=motion.record, level=motion.level, period=2 * math.pi / omega,
        time=u[:n, 0], u=u[:n, 1], force=-frc[:n, 0] if frc.size else np.zeros(n),
        completed=end_time >= t_max - 0.5 * dt, end_time=end_time, duration=t_max, collapsed=collapsed,
    )


CYCLIC_PROTOCOL = [1, 0, -1, 0, 2, 0, -2, 0, 5, 0, -5, 0, 8, 0, -8, 0, 10, 0, -10, 0, 15, 0, -15, 0, 20, 0, -20, 0,
                   25, 0, -25, 0, 30, 0, -30, 0, 35, 0, -35, 0]


@dataclass
class SDOFPushover:
    u: np.ndarray           # SDOF displacement (mm)
    force: np.ndarray       # base shear = - support reaction (N)
    completed: bool


def run_sdof_pushover(model: SDOFModel, targets: list[float] | None = None, step: float = 0.05,
                      tolerance: float = 1e-4, max_iterations: int = 1000, fallback_tolerance: float = 1e-3,
                      fallback_max_iterations: int = 5000) -> SDOFPushover:
    """Displacement-controlled pushover of the SDOF (springs without caps), as ``Pushover_SDOF/<tag>_SDOF.py``:
    monotonic to ``targets = [50]`` mm by default, or through a list of target displacements (e.g.
    ``CYCLIC_PROTOCOL``), in steps of ``step`` with the script's fallbacks."""
    targets = [50.0] if targets is None else list(targets)
    _build(model)
    op.eigen(1)
    with tempfile.TemporaryDirectory() as tmp:
        f_r, f_u = Path(tmp) / "r.out", Path(tmp) / "u.out"
        op.recorder("Node", "-file", str(f_r), "-node", 1, "-dof", 1, "reaction")
        op.recorder("Node", "-file", str(f_u), "-node", 2, "-dof", 1, "disp")
        op.numberer("RCM")
        op.system("FullGeneral")
        op.constraints("Transformation")
        op.test("EnergyIncr", tolerance, max_iterations)
        op.algorithm("Newton")
        op.analysis("Static")
        op.timeSeries("Linear", 2)
        op.pattern("Plain", 100, 2)
        op.load(2, 1, 0, 0)
        op.reactions()

        ok = 0
        for j, target in enumerate(targets):
            if ok != 0:
                break
            disp = target - (targets[j - 1] if j > 0 else 0)
            dstep = step if disp > 0 else -step
            op.integrator("DisplacementControl", 2, 1, dstep)
            for _ in range(int(disp / dstep)):
                ok = op.analyze(1)
                for algorithm, n_sub in ((("Newton", "-initial"), 2), (("Broyden", 50), 2), (("NewtonLineSearch",), 10)):
                    if ok == 0:
                        break
                    op.test("EnergyIncr", fallback_tolerance, fallback_max_iterations)
                    op.algorithm(*algorithm)
                    op.integrator("DisplacementControl", 2, 1, dstep / n_sub)
                    ok = op.analyze(n_sub)
                    op.integrator("DisplacementControl", 2, 1, dstep)
                    op.test("EnergyIncr", tolerance, max_iterations)
                    op.algorithm("Newton")
                if ok != 0:
                    break
        op.wipe()
        r = np.loadtxt(f_r, ndmin=1)
        u = np.loadtxt(f_u, ndmin=1)
    n = min(len(r), len(u))
    return SDOFPushover(u=u[:n], force=-r[:n], completed=ok == 0)


def support_demands(model: SDOFModel, peak_u: float) -> list[dict]:
    """Peak support displacement Gamma * phi * u_SDOF and its ratio to the trapeze envelope points."""
    if model.gamma is None:
        return []
    rows = []
    for label, phi, spring in zip(model.support_labels, model.support_phi, model.springs):
        d = model.gamma * phi * peak_u
        base = [v * model.gamma * phi for v in spring.pos_disp]      # unscaled trapeze deformations
        rows.append({"support": label, "phi": phi, "peak_displacement": d,
                     **{f"ratio_to_ePd{i + 1}": d / base[i] for i in range(4)}})
    return rows
