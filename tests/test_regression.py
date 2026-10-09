"""The engine against the scripts in ``code implementation for paper/`` and their stored results."""

import numpy as np
import pytest

from piperom.inputs import AnalysisSettings, load_system
from piperom.sdof import derive_sdof
from validation.legacy import (ARCHETYPES, EQUIV_STATIC_EXAMPLE, archetype, compare_inputs, compare_pushover,
                               extract_support_displacements, procedure_sdof, procedure_settings, procedure_static,
                               script_sdof)


@pytest.mark.parametrize("tag", ARCHETYPES)
def test_archetype_inputs_match_the_pushover_drivers(tag):
    assert compare_inputs(tag) == []


@pytest.mark.parametrize("tag", ARCHETYPES)
def test_pushover_reproduces_the_scripts(tag):
    """Every value of Pushover2D/pushover_results_<tag>.txt (written with 3 decimals) is reproduced."""
    cmp = compare_pushover(tag)
    assert cmp["n_diff"] == 0, f"{tag}: {cmp['n_diff']} of {cmp.get('n_cells')} values differ"
    assert cmp["max_abs_diff"] <= 0.0005 + 1e-9


@pytest.mark.parametrize("tag", ARCHETYPES)
def test_sdof_reproduces_sdof_from_2d(tag):
    """derive_sdof = Pushover_SDOF/sdof_from_2d.py (which reads the pushover results, 3 decimals)."""
    ref = script_sdof(tag)
    p = derive_sdof(archetype(tag), AnalysisSettings.from_dict())
    assert (p.delta_c, p.gamma, p.effective_mass) == (ref["dc"], ref["Gamma"], ref["mass"])
    assert [(s.phi, s.n_trapezes) for s in p.supports if s.kind == "transverse"] == [tuple(v) for v in ref["trans"]]
    assert [(s.phi, s.n_trapezes) for s in p.supports if s.kind == "longitudinal"] == [tuple(v) for v in ref["long"]]


def test_sdof_reproduces_the_proposed_procedure_sdof():
    """derive_sdof = the SDOF NLTHA_SDOF.py analyses (Results/NLTHA/SDOF_nT4_nL3/sdof_params.json)."""
    ref = procedure_sdof()
    p = derive_sdof(load_system(EQUIV_STATIC_EXAMPLE), procedure_settings(), 12.0)
    assert (p.delta_c, p.gamma, p.effective_mass) == (ref["dc"], ref["Gamma"], ref["mass"])
    assert [[s.phi, s.n_trapezes] for s in p.supports if s.kind == "transverse"] == ref["trans"]
    assert [[s.phi, s.n_trapezes] for s in p.supports if s.kind == "longitudinal"] == ref["long"]


def test_sdof_reproduces_equivalent_static_procedure():
    """derive_sdof at full precision = code_proposed_procedure/equivalent_static.py (stored at full precision)."""
    ref = procedure_static()
    rs = load_system(EQUIV_STATIC_EXAMPLE).resolve()
    params = derive_sdof(rs, procedure_settings(round_decimals=None), delta_c=float(ref["dc"]))

    assert np.array_equal(rs.brace_mask, ref["stiff_mask"])
    assert params.n_transverse == int(ref["nT"])
    assert params.n_longitudinal_trapezes == int(ref["nL"])
    assert params.gamma == pytest.approx(float(ref["Gamma"]), rel=1e-12)
    assert params.effective_mass == pytest.approx(float(ref["M_eff"]), rel=1e-12)
    assert params.u_sdof == pytest.approx(float(ref["u_sdof"]), rel=1e-12)
    np.testing.assert_allclose(params.d_norm, ref["d_norm"], rtol=1e-12)
    disp_shape = extract_support_displacements(ref["d_norm"], ref["stiff_mask"], int(ref["nOrth"]))
    np.testing.assert_allclose(params.support_shape, disp_shape, rtol=1e-12)

    # Scaled Pinching4 envelopes as built in NLTHA_SDOF.py: forces x n, deformations / (Gamma phi)
    gamma = float(ref["Gamma"])
    for s, phi in zip(params.supports, disp_shape):
        base = rs.transverse if s.kind == "transverse" else rs.longitudinal
        n = 1 if s.kind == "transverse" else int(ref["nL"])
        np.testing.assert_allclose(s.spring.pos_disp, [v / (gamma * phi) for v in base.pos_disp], rtol=1e-12)
        np.testing.assert_allclose(s.spring.neg_disp, [v / (gamma * phi) for v in base.neg_disp], rtol=1e-12)
        np.testing.assert_allclose(s.spring.pos_force, [n * v for v in base.pos_force], rtol=1e-12)
        np.testing.assert_allclose(s.spring.neg_force, [n * v for v in base.neg_force], rtol=1e-12)
