"""Helpers comparing the engine with the scripts and stored results in ``code implementation for paper/``.

The reference is the scripts' current state (commit c13ed81 and the corrections of 2026-10-09): the pushover drivers
``Pushover2D/CS_Lumped_Iter_<tag>_cont_PO.py`` and their ``pushover_results_<tag>.txt``, the SDOFs of
``Pushover_SDOF/sdof_from_2d.py`` and their time histories ``Pushover_SDOF/Results<tag>_SDOF/``, and the
3D models ``3D_models/<M>_biron.py`` with ``Results/<M>_biron_DispX/Y.txt``.

``code_proposed_procedure/`` applies the same model to one system at one Delta_c
(``equivalent_static.py``, results in ``Results/EquivStatic/``) and runs its SDOF under the 150 floor motions
(``NLTHA_SDOF.py``, results in ``Results/NLTHA/SDOF_nT4_nL3/``).
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from piperom.inputs import AnalysisSettings, PipingSystem, load_system  # noqa: E402
from piperom.pushover import PushoverResult, run_pushover  # noqa: E402

PAPER = REPO / "code implementation for paper"   # original scripts and results
ARCHETYPE_DIR = REPO / "inputs" / "archetypes"
ARCHETYPES = sorted(p.stem for p in ARCHETYPE_DIR.glob("*.yaml"))
EQUIV_STATIC_EXAMPLE = REPO / "inputs" / "examples" / "equivalent_static_example.yaml"
PROCEDURE = PAPER / "code_proposed_procedure" / "Results"
PROCEDURE_TAG = "equiv_static_nT4_nL3"


def archetype(tag: str) -> PipingSystem:
    return load_system(ARCHETYPE_DIR / f"{tag}.yaml")


def procedure_settings(round_decimals: int | None = 3) -> AnalysisSettings:
    """Settings of ``code_proposed_procedure/``: the SDOF at exactly Delta_c (its only static step); with
    ``round_decimals=None`` the SDOF keeps full precision (as stored in the ``.npz``)."""
    return AnalysisSettings.from_dict({"sdof": {"on_pushover_grid": False, "round_decimals": round_decimals}})


def procedure_static() -> dict:
    """Stored output of ``code_proposed_procedure/equivalent_static.py`` (full precision)."""
    data = np.load(PROCEDURE / "EquivStatic" / f"{PROCEDURE_TAG}.npz")
    return {k: data[k] for k in data.files}


def procedure_sdof() -> dict:
    """``Results/NLTHA/SDOF_nT4_nL3/sdof_params.json``: the SDOF ``NLTHA_SDOF.py`` analysed (3 decimals)."""
    return json.loads((PROCEDURE / "NLTHA" / "SDOF_nT4_nL3" / "sdof_params.json").read_text())


def procedure_peaks() -> list[float]:
    """Peak |u| of the SDOF under each S4_150 record (``Disp_GM<k>_Floor4.txt``: min, max, abs max)."""
    d = PROCEDURE / "NLTHA" / "SDOF_nT4_nL3"
    return [float(np.loadtxt(d / f"Disp_GM{k}_Floor4.txt")[2]) for k in range(1, 151)]


def legacy_driver(tag: str) -> Path:
    return PAPER / "Pushover2D" / f"CS_Lumped_Iter_{tag}_cont_PO.py"


class _Stop(Exception):
    pass


def driver_inputs(tag: str) -> dict:
    """Keyword arguments the pushover driver of ``tag`` passes to ``iterate_shape_from_static``."""
    sys.path.insert(0, str(PAPER / "Pushover2D"))
    captured: dict = {}

    def stub(**kw):
        captured.update(kw)
        raise _Stop

    ns: dict = {"__name__": "__capture__", "__file__": str(legacy_driver(tag))}
    exec(compile("from Functions import *\n", "<capture>", "exec"), ns)
    ns["iterate_shape_from_static"] = stub
    src = legacy_driver(tag).read_text(encoding="utf-8").replace("from Functions import*", "")
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(src, str(legacy_driver(tag)), "exec"), ns)
    except _Stop:
        pass
    return captured


def compare_inputs(tag: str) -> list[str]:
    """Differences between ``inputs/archetypes/<tag>.yaml`` and the driver's inputs (empty if none)."""
    kw = driver_inputs(tag)
    rs = archetype(tag).resolve()
    p = rs.pipe
    checks = {
        "length": (rs.length, kw["Lpipe"]),
        "n_pipes": (rs.n_pipes, kw["npipes"]),
        "n_mains_left/right": ((rs.n_mains_left, rs.n_mains_right), (kw["n_mains_left"], kw["n_mains_right"])),
        "x_center": (rs.x_center, kw["x_center"]),
        "pipe": ((p.outer_diameter, p.inner_diameter, p.elastic_modulus, p.shear_modulus, p.density),
                 (kw["Dext"], kw["Dint"], kw["E"], kw["G"], kw["rho"])),
        "hangers": (rs.hanger_x.tolist(),
                    [float(x) for x in kw["x_hangers_user"]] if kw.get("x_hangers_user") is not None
                    else [float(x) for x in kw["x0_soft"] + kw["soft_spacing"] * np.arange(len(rs.hanger_x))]),
        "brace mask": (rs.brace_mask.tolist(), np.asarray(kw["stiff_mask"]).tolist()),
        "branches": ([rs.branch_x.tolist(), rs.branch_length.tolist(), rs.branch_n_pipes.tolist(),
                      rs.branch_n_braces.tolist()],
                     [np.asarray(kw[k], float).tolist() for k in
                      ("x_ortho_user", "L_ortho_user", "n_ortho_pipes_user", "n_ortho_springs_user")]),
        "branch participation": (rs.branch_participation, kw.get("alpha", 1.0)),
    }
    return [f"{k}: yaml {a} vs driver {b}" for k, (a, b) in checks.items() if a != b]


def legacy_table(result: PushoverResult, tag: str) -> np.ndarray:
    """Engine results laid out like ``pushover_results_<tag>.txt`` (incl. the driver's DOF sorting)."""
    sort = "d_norm[order]" in legacy_driver(tag).read_text(encoding="utf-8")
    order = np.argsort(np.array([d["x"] for d in result.dofs])) if sort else slice(None)
    return np.array([[s.delta_c, s.gamma, s.effective_mass, s.base_shear, s.mass_ratio, s.u_sdof,
                      *s.d_norm[order], *s.d_scaled[order], *s.f_push[order]] for s in result.steps])


def compare_pushover(tag: str) -> dict:
    """Engine (default settings) vs ``Pushover2D/pushover_results_<tag>.txt`` (3 decimals)."""
    result = run_pushover(archetype(tag), AnalysisSettings.from_dict())
    new = legacy_table(result, tag)
    ref = np.loadtxt(PAPER / "Pushover2D" / f"pushover_results_{tag}.txt")
    if new.shape != ref.shape:
        return {"tag": tag, "shape_new": new.shape, "shape_ref": ref.shape, "n_diff": -1}
    formatted = np.vectorize(lambda v: float(f"{v:.3f}"))(new)
    return {"tag": tag, "n_cells": ref.size, "n_diff": int((formatted != ref).sum()),
            "max_abs_diff": float(np.abs(new - ref).max()), "result": result}


def script_sdof(tag: str, dc_target: float = 12.0) -> dict:
    """``load_sdof_params`` of ``Pushover_SDOF/sdof_from_2d.py``: the SDOF the NLTHA scripts analyse."""
    sys.path.insert(0, str(PAPER / "Pushover_SDOF"))
    from sdof_from_2d import load_sdof_params
    with contextlib.redirect_stdout(io.StringIO()):
        return load_sdof_params(tag, dc_target)


def script_sdof_model(tag: str | None = None, sd: dict | None = None):
    """The SDOF of ``Pushover_SDOF/<tag>_SDOF_NLTHA.py`` (or the one given as ``sd``, in the format of
    ``load_sdof_params``) as an engine model, with the values the script reads."""
    from piperom.timehistory import SDOFModel
    from piperom.trapeze import load_trapeze
    sd = script_sdof(tag) if sd is None else sd
    tag = tag or "sdof"
    t, l = load_trapeze("default", "transverse"), load_trapeze("default", "longitudinal")
    g = sd["Gamma"]
    springs = [l.scaled(n, g * phi) for phi, n in sd["long"]] + [t.scaled(n, g * phi) for phi, n in sd["trans"]]
    return SDOFModel(name=tag, mass=sd["mass"], springs=springs, gamma=g,
                     support_phi=[phi for phi, _ in sd["long"] + sd["trans"]],
                     support_kinds=["longitudinal"] * len(sd["long"]) + ["transverse"] * len(sd["trans"]))


def script_sdof_peaks(tag: str) -> dict[tuple[int, str], tuple[float, bool, bool]]:
    """``Pushover_SDOF/Results<tag>_SDOF/<tag>_NLTHA_peak_displacements.csv``:
    (IM, record) -> (peak |u| in mm, completed, collapsed)."""
    import csv
    path = PAPER / "Pushover_SDOF" / f"Results{tag}_SDOF" / f"{tag}_NLTHA_peak_displacements.csv"
    with open(path, newline="") as f:
        return {(int(r["IM"]), r["record"]): (float(r["peak_disp_mm"]), r["completed"] == "1", r["collapsed"] == "1")
                for r in csv.DictReader(f)}


def script_3d_peak(model: str, direction: str, record_index: int, level: int) -> float:
    """``Results/<model>_biron_Disp<X|Y>.txt``: largest peak over the braced nodes (100 = collapse, 0 = not
    completed); row = index of the record applied in that direction, column = IM."""
    return float(np.loadtxt(PAPER / "Results" / f"{model}_biron_Disp{direction.upper()}.txt")[record_index, level - 1])


def stored_3d_peak(peak: float, completed: bool, collapsed: bool) -> float:
    """A 3D result as the scripts store it (``3D_models/*_biron.py``, 2 decimals)."""
    return 100.0 if collapsed else (round(peak, 2) if completed else 0.0)


def extract_support_displacements(d_dof, stiff_mask, n_orth):
    """Copy of the function in ``code_proposed_procedure/NLTHA_SDOF.py``."""
    d_dof = np.asarray(d_dof)
    stiff_mask = np.asarray(stiff_mask, dtype=bool)
    d_stiff = d_dof[:len(stiff_mask)][stiff_mask]
    return np.concatenate([d_stiff, d_dof[-n_orth:]]) if n_orth > 0 else d_stiff


# ----------------------------------------------------------------------------- capturing OpenSees models

class StopCapture(Exception):
    pass


class OpenSeesRecorder:
    """Stand-in for ``openseespy.opensees`` that records calls instead of executing them."""

    def __init__(self, stop_at: set[str]):
        self.calls: list[tuple[str, tuple]] = []
        self._stop_at = stop_at

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)

        def call(*args):
            if name in self._stop_at:
                raise StopCapture(name)
            self.calls.append((name, args))
            return [] if name == "eigen" else 0
        return call


def capture_opensees_calls(script: Path, stop_at: set[str], patch=None) -> tuple[list, dict]:
    """Execute ``script`` with a recording ``op`` until one of ``stop_at`` is called.

    Returns the recorded calls and the script's global variables at that point. The script runs in a
    temporary working directory, so folders it creates don't end up in the repository.
    """
    import os
    src = script.read_text(encoding="utf-8")
    src = src.replace("import openseespy.opensees as op", "").replace("import vfo.vfo as vfo", "")
    if patch:
        src = patch(src)
    rec = OpenSeesRecorder(stop_at)
    ns = {"__name__": "__capture__", "__file__": str(Path(script).resolve()), "op": rec}
    cwd, argv = os.getcwd(), sys.argv
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)
        sys.argv = [str(script)]          # scripts reading command-line options get their defaults
        try:
            exec(compile(src, str(script), "exec"), ns)
        except StopCapture:
            pass
        finally:
            os.chdir(cwd)
            sys.argv = argv
    return rec.calls, ns
