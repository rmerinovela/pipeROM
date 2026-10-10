"""Time-history engines against the scripts' stored results (need the floor motions in motions/)."""

import pytest

from piperom.inputs import InputError
from piperom.motions import MOTIONS_DIR, load_motion_sets, motion_set, select_runs
from piperom.timehistory import TimeHistorySettings, run_sdof_time_history
from piperom.trapeze import Pinching4, load_trapeze
from validation.legacy import (procedure_peaks, script_3d_peak, script_sdof_model,
                               script_sdof_peaks, stored_3d_peak)

HAVE_MOTIONS = (MOTIONS_DIR / "floor_motions" / "ResultsS4" / "FloorAcc_IM10_120111.txt").exists()
needs_motions = pytest.mark.skipif(not HAVE_MOTIONS, reason="floor motions not available in motions/")


def test_pinching4_opensees_args_roundtrip():
    p = load_trapeze("default", "longitudinal")
    again = Pinching4.from_opensees_args(p.opensees_args(), name=p.name)
    assert again == p


def test_motion_catalogue_and_selection():
    sets = load_motion_sets()
    assert {"S4_IM", "S4_150"} <= set(sets)
    ms = sets["S4_IM"]
    assert len(ms.records) == 44 and len(ms.pairs()) == 22 and ms.levels == list(range(1, 11))
    assert len(select_runs(ms, [1, 2], None)) == 88
    assert select_runs(ms, [3], ["120111"]) == [("120111", 3)]
    assert len(select_runs(sets["S4_150"])) == 150
    with pytest.raises(InputError, match="not available"):
        select_runs(ms, [13], None)
    with pytest.raises(InputError, match="not in motion set"):
        select_runs(ms, [1], ["999"])


@needs_motions
@pytest.mark.parametrize("tag, cases", [
    ("M01x", [(1, "120111"), (10, "120112"), (5, "120411")]),
    ("M03x", [(10, "120111"), (7, "120112"), (3, "120411")]),     # M03x collapses at the strong IMs
    ("M62y", [(10, "120111"), (10, "120411"), (4, "120112")]),
])
def test_sdof_time_history_reproduces_the_scripts(tag, cases):
    """Pushover_SDOF/<tag>_SDOF_NLTHA.py: peak |u| (4 decimals), completed and collapsed flags, with the
    engine's own SDOF (default settings)."""
    from piperom.inputs import AnalysisSettings
    from piperom.sdof import derive_sdof
    from piperom.timehistory import SDOFModel
    from validation.legacy import archetype
    ref = script_sdof_peaks(tag)
    model = SDOFModel.from_parameters(derive_sdof(archetype(tag), AnalysisSettings.from_dict()))
    ms, s = motion_set("S4_IM"), TimeHistorySettings.from_dict()
    for level, record in cases:
        r = run_sdof_time_history(model, ms.load(record, level), s)
        peak, completed, collapsed = ref[(level, record)]
        assert f"{r.peak_u:.4f}" == f"{peak:.4f}" and (r.completed, r.collapsed) == (completed, collapsed)


@pytest.mark.parametrize("tag", ["M01x", "M29x", "M62y"])
def test_sdof_pushover_reproduces_the_scripts(tag):
    """Pushover_SDOF/<tag>_SDOF.py: DispC.out and VbaseC.out (support reaction), 6 significant digits.
    (M61y's stored files predate its current pushover results.)"""
    import numpy as np
    from piperom.timehistory import run_sdof_pushover
    from validation.legacy import PAPER
    d = PAPER / "Pushover_SDOF" / f"Results{tag}_SDOF"
    po = run_sdof_pushover(script_sdof_model(tag))
    assert po.completed
    for new, ref in ((po.u, np.loadtxt(d / "DispC.out")), (-po.force, np.loadtxt(d / "VbaseC.out"))):
        assert [f"{v:.6g}" for v in new] == [f"{v:.6g}" for v in ref]


@needs_motions
def test_proposed_procedure_time_history_reproduces_the_script():
    """code_proposed_procedure/NLTHA_SDOF.py under S4_150 (first records), with the engine's own SDOF."""
    from piperom.sdof import derive_sdof
    from piperom.timehistory import SDOFModel
    from validation.legacy import EQUIV_STATIC_EXAMPLE, procedure_settings
    from piperom.inputs import load_system
    model = SDOFModel.from_parameters(derive_sdof(load_system(EQUIV_STATIC_EXAMPLE), procedure_settings(), 12.0))
    ms, s = motion_set("S4_150"), TimeHistorySettings.from_dict()
    peaks = procedure_peaks()
    for k in range(3):
        r = run_sdof_time_history(model, ms.load(ms.records[k], None, 4), s)
        assert f"{r.peak_u:.4f}" == f"{peaks[k]:.4f}"


@needs_motions
def test_3d_verification_reproduces_the_scripts():
    """Results/M01_biron_DispX/Y.txt: largest peak over the braced nodes, 2 decimals (row = record index)."""
    from piperom.verification3d import Verification3DSettings, braced_peaks, load_model3d, run_3d
    ms, m = motion_set("S4_IM"), load_model3d("M01")
    rx, ry = ms.pairs()[0]
    r = run_3d(m, ms.load(rx, 10), ms.load(ry, 10), Verification3DSettings.from_dict())
    for d, peak in braced_peaks(m, r).items():
        rec = rx if d == "x" else ry
        ref = script_3d_peak("M01", d, ms.records.index(rec), 10)
        assert f"{stored_3d_peak(peak, r.completed, r.collapsed):.2f}" == f"{ref:.2f}"
