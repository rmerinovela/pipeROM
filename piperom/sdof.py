"""Equivalent SDOF parameters at a chosen target displacement Delta_c.

Reproduces ``code implementation for paper/Pushover_SDOF/sdof_from_2d.py`` (``load_sdof_params``):

* the equivalent static procedure is run at Delta_c; by default (``sdof.on_pushover_grid``) Delta_c is
  the pushover step closest to the requested value, as the scripts read the SDOF from the pushover results,
  and the values are rounded to 3 decimals as written there (``sdof.round_decimals``);
* Gamma and the effective mass come from the converged shape;
* each braced hanger (transverse) and each branch (longitudinal) is one support, with
  phi = displaced shape at the support, normalised to the reference (last) branch DOF;
* each support becomes a Pinching4 spring of the SDOF with envelope deformations divided by
  Gamma * phi and envelope forces multiplied by the number of trapezes it represents: the number of
  mains lumped on its side of x_center for a braced hanger, ``alpha * n_braces`` for a branch. The SDOF
  base shear then matches the base shear of the static model.

The SDOF itself (mass = effective mass, parallel springs) is analysed in ``timehistory``.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .inputs import AnalysisSettings, InputError, PipingSystem, ResolvedSystem, dump_yaml
from .pushover import PushoverStep, run_step
from .trapeze import Pinching4

MIN_PHI = 1e-3      # a spring needs a positive shape value at its support


@dataclass
class SDOFSupport:
    kind: str           # "transverse" (braced hanger) or "longitudinal" (branch)
    dof: int
    x: float
    branch: int | None
    phi: float
    n_trapezes: float   # trapezes represented by the spring (force scale)
    spring: Pinching4   # scaled Pinching4 of the SDOF spring


@dataclass
class SDOFParameters:
    name: str
    delta_c: float
    gamma: float
    effective_mass: float
    total_mass: float
    mass_ratio: float
    u_sdof: float
    base_shear: float
    converged: bool
    iterations: int
    supports: list[SDOFSupport]
    d_norm: np.ndarray

    @property
    def n_transverse(self) -> int:
        return sum(s.kind == "transverse" for s in self.supports)

    @property
    def n_longitudinal_trapezes(self) -> float:
        return sum(s.n_trapezes for s in self.supports if s.kind == "longitudinal")

    @property
    def support_shape(self) -> np.ndarray:
        """phi at transverse supports then branches (``DispShape`` of the original SDOF scripts)."""
        return np.array([s.phi for s in self.supports])

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "delta_c": self.delta_c,
            "gamma": self.gamma,
            "effective_mass": self.effective_mass,
            "total_mass": self.total_mass,
            "mass_ratio": self.mass_ratio,
            "u_sdof": self.u_sdof,
            "base_shear": self.base_shear,
            "converged": self.converged,
            "iterations": self.iterations,
            "n_transverse_supports": self.n_transverse,
            "n_longitudinal_trapezes": self.n_longitudinal_trapezes,
            "supports": [
                {"kind": s.kind, "dof": s.dof, "x": s.x, "branch": s.branch, "phi": float(s.phi),
                 "n_trapezes": s.n_trapezes, "pinching4": s.spring.parameters()}
                for s in self.supports
            ],
        }

    def to_yaml(self) -> str:
        return dump_yaml(self.to_dict())

    def springs_csv(self) -> str:
        """One row per SDOF spring with its scaled envelope (forces in N, deformations in mm of SDOF)."""
        buf = io.StringIO()
        cols = (["support", "kind", "x", "branch", "phi", "n_trapezes"]
                + [f"ePf{i}" for i in range(1, 5)] + [f"ePd{i}" for i in range(1, 5)]
                + [f"eNf{i}" for i in range(1, 5)] + [f"eNd{i}" for i in range(1, 5)])
        w = csv.writer(buf, lineterminator="\n")
        w.writerow(cols)
        for k, s in enumerate(self.supports):
            p = s.spring.parameters()
            w.writerow([k, s.kind, repr(s.x), "" if s.branch is None else s.branch, repr(float(s.phi)),
                        repr(float(s.n_trapezes))] + [repr(float(p[c])) for c in cols[6:]])
        return buf.getvalue()

    def write(self, out_dir: str | Path) -> list[Path]:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        files = [out / "sdof_parameters.yaml", out / "sdof_springs.csv"]
        files[0].write_text(self.to_yaml())
        files[1].write_text(self.springs_csv())
        return files


def _rounded(v, decimals: int | None):
    """``v`` as written with ``decimals`` decimals and read back (``np.savetxt(fmt="%.3f")`` + ``np.loadtxt``)."""
    if decimals is None:
        return v
    if isinstance(v, np.ndarray):
        return np.array([float(f"{x:.{decimals}f}") for x in v])
    return float(f"{v:.{decimals}f}")


def sdof_from_step(rs: ResolvedSystem, step: PushoverStep, converged: bool | None = None,
                   iterations: int | None = None, decimals: int | None = None) -> SDOFParameters:
    """SDOF from a pushover step. With ``decimals``, the values of the step's row of the pushover results
    (Delta_c, Gamma, M_eff, V_b, mass ratio, u_SDOF, shape) are rounded as the scripts read them."""
    r = lambda v: _rounded(v, decimals)  # noqa: E731
    gamma = r(step.gamma)
    d_norm = r(step.d_norm)
    supports: list[SDOFSupport] = []
    for i in np.where(rs.brace_mask == 1)[0]:
        x = float(rs.hanger_x[i])
        phi = float(d_norm[i])
        n = rs.n_mains_at(x)
        supports.append(SDOFSupport("transverse", int(i), x, None, phi, n,
                                    rs.transverse.scaled(n, gamma * phi)))
    for j in range(rs.n_branches):
        i = rs.n_hangers + j
        phi = float(d_norm[i])
        n = rs.branch_participation * int(rs.branch_n_braces[j])
        supports.append(SDOFSupport("longitudinal", i, float(rs.branch_x[j]), j, phi, n,
                                    rs.longitudinal.scaled(n, gamma * phi)))
    for s in supports:
        if s.phi <= MIN_PHI:
            raise InputError(f"{rs.name}: non-positive shape value {s.phi:.3g} at the {s.kind} support at "
                             f"x = {s.x:g} (Delta_c = {step.delta_c:.2f} mm); choose another Delta_c")
    return SDOFParameters(
        name=rs.name, delta_c=r(step.delta_c), gamma=gamma, effective_mass=r(step.effective_mass),
        total_mass=step.total_mass, mass_ratio=r(step.mass_ratio), u_sdof=r(step.u_sdof),
        base_shear=r(step.base_shear),
        converged=step.converged if converged is None else converged,
        iterations=step.iterations if iterations is None else iterations,
        supports=supports, d_norm=d_norm,
    )


def sdof_delta_c(settings: AnalysisSettings, delta_c: float | None = None) -> float:
    """Delta_c at which the SDOF is derived: the requested value (default ``sdof.delta_c``), or the pushover
    step closest to it when ``sdof.on_pushover_grid`` is set."""
    dc = settings.sdof_delta_c if delta_c is None else float(delta_c)
    if settings.sdof_on_pushover_grid:
        grid = settings.delta_c()
        dc = float(grid[np.argmin(np.abs(grid - dc))])
    return dc


def derive_sdof(system: PipingSystem | ResolvedSystem, settings: AnalysisSettings,
                delta_c: float | None = None) -> SDOFParameters:
    """Run the equivalent static procedure at ``sdof_delta_c(settings, delta_c)``."""
    rs = system.resolve() if isinstance(system, PipingSystem) else system
    step, _ = run_step(rs, settings, sdof_delta_c(settings, delta_c))
    return sdof_from_step(rs, step, decimals=settings.sdof_round_decimals)
