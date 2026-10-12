"""First-yield step of the pseudo-pushover (the displaced shape shown by default in the app)."""

import numpy as np

from piperom.inputs import AnalysisSettings
from piperom.pushover import first_yield_step, run_pushover
from validation.legacy import archetype


def restraint_ratio(rs, step):
    """Largest restraint displacement over its first-yield displacement."""
    braced = np.where(rs.brace_mask == 1)[0]
    D = np.abs(step.d_scaled)
    return max(np.max(D[braced]) / rs.transverse_trilinear.d1,
               np.max(D[rs.n_hangers:]) / rs.longitudinal_trilinear.d1)


def test_first_yield_step_is_the_first_step_with_a_yielded_restraint():
    rs = archetype("M02y").resolve()
    settings = AnalysisSettings.from_dict({"pushover": {"delta_c_values": [6, 9, 11, 12, 13, 16]}})
    res = run_pushover(rs, settings)
    k = first_yield_step(res, rs)
    ratios = [restraint_ratio(rs, s) for s in res.steps]
    assert k is not None and ratios[k] >= 1.0
    assert all(r < 1.0 for r in ratios[:k])


def test_no_first_yield_step_in_the_elastic_range():
    rs = archetype("M02y").resolve()
    res = run_pushover(rs, AnalysisSettings.from_dict({"pushover": {"delta_c_values": [2, 4]}}))
    assert first_yield_step(res, rs) is None
