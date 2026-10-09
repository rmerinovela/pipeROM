"""Equilibrium of the equivalent static load pattern and the base shear (docs/legacy_issues.md, A2)."""

import numpy as np
import pytest

from piperom.static_model import SolverSettings, solve_static_step
from validation.legacy import ARCHETYPES, archetype


def static_step(tag, split, delta=20.0):
    rs = archetype(tag).resolve()
    d = np.ones(rs.n_dof)
    step = solve_static_step(rs, d, delta, SolverSettings("NormDispIncr", 1e-8, 50, split))
    braced = np.where(rs.brace_mask == 1)[0]
    brace_f = step.brace_stiffness * delta * d[braced]
    branch_f = step.branch_stiffness * delta * d[rs.n_hangers:]
    return rs, step, brace_f, branch_f


@pytest.mark.parametrize("tag", ARCHETYPES)
def test_consistent_split_applies_the_full_spring_forces(tag):
    """The load pattern of the analysed main equals its support forces at the assumed shape."""
    _, step, brace_f, branch_f = static_step(tag, "consistent")
    applied = np.sum(step.loads) + np.sum(step.branch_stay_load)
    assert applied == pytest.approx(np.sum(brace_f) + np.sum(branch_f), rel=1e-12)


@pytest.mark.parametrize("tag", ARCHETYPES)
def test_base_shear_is_the_support_reactions_of_all_mains(tag):
    rs, step, brace_f, branch_f = static_step(tag, "consistent")
    n_brace = [rs.n_mains_at(x) for x in rs.hanger_x[rs.brace_mask == 1]]
    n_branch = [rs.n_mains_at(x) for x in rs.branch_x]
    assert step.base_shear == pytest.approx(np.dot(n_brace, brace_f) + np.dot(n_branch, branch_f), rel=1e-12)


def test_legacy_split_loses_part_of_the_branch_force():
    _, step, brace_f, branch_f = static_step("M01y", "legacy")
    assert np.sum(step.loads) + np.sum(step.branch_stay_load) < np.sum(brace_f) + np.sum(branch_f)


def test_lumped_mains():
    """M29x lumps 6 mains left of x_center and 9 right of it: a massless node at x_center, bending
    stiffness and modal mass scaled per side, branch braces shared by the mains."""
    rs, step, _, _ = static_step("M29x", "consistent")
    assert (rs.n_mains_left, rs.n_mains_right, rs.x_center) == (6, 9, 5000.0)
    assert rs.n_mains_at(500) == 6 and rs.n_mains_at(9500) == 9 and rs.n_mains_at(5000) == 7.5
    n = np.array([rs.n_mains_at(x) for x in step.dof_x])
    assert step.total_mass == pytest.approx(np.sum(n * step.dof_mass) + np.sum(step.branch_mass), rel=1e-12)
    # a single main gives the same tributary masses: the x_center node carries no mass
    one = archetype("M29x")
    one.n_mains_left = one.n_mains_right = None
    one.n_mains = 1
    _, step1, _, _ = static_step_from(one.resolve())
    np.testing.assert_allclose(step1.dof_mass, step.dof_mass, rtol=1e-12)


def static_step_from(rs, delta=20.0):
    d = np.ones(rs.n_dof)
    step = solve_static_step(rs, d, delta, SolverSettings("NormDispIncr", 1e-8, 50, "consistent"))
    return rs, step, None, None
