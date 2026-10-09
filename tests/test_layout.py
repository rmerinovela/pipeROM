"""Plan layouts and their reduction to the equivalent system of each direction (piperom.layout)."""

import copy

import numpy as np
import pytest

from piperom.inputs import InputError
from piperom.layout import layout_from_dict, load_layout
from validation.legacy import REPO, archetype

LAYOUTS = REPO / "inputs" / "layouts"
# Directions whose 2D model in the paper is consistent with its 3D model (docs/legacy_issues.md, E5-E6)
CONSISTENT = ["M01x", "M01y", "M02x", "M02y", "M03x", "M03y", "M29x", "M29y", "M30x", "M30y", "M31x", "M31y",
              "M61x", "M61y", "M62x", "M62y", "M63x", "M63y"]


def resolved_arrays(system):
    rs = system.resolve()
    return (rs.length, rs.n_pipes, rs.n_mains_left, rs.n_mains_right, rs.x_center, rs.hanger_x.tolist(),
            rs.brace_mask.tolist(), rs.branch_x.tolist(), rs.branch_length.tolist(), rs.branch_n_pipes.tolist(),
            rs.branch_n_braces.tolist())


@pytest.mark.parametrize("tag", CONSISTENT)
def test_layout_reduces_to_the_archetype_files(tag):
    system = load_layout(LAYOUTS / f"{tag[:3]}.yaml").systems()[tag[3]]
    assert resolved_arrays(system) == resolved_arrays(archetype(tag))


def base():
    return load_layout(LAYOUTS / "M29.yaml").raw


def test_branches_must_have_the_same_length():
    d = copy.deepcopy(base())
    d["branches"][0]["length"] = 6000
    d["branches"][0]["hangers"] = {"positions": [1500, 4500]}
    with pytest.raises(InputError, match="same length"):
        layout_from_dict(d)


def test_branches_on_a_side_must_be_identical():
    d = copy.deepcopy(base())
    d["branches"][0]["braces"]["transverse"] = {"positions": [4500]}
    with pytest.raises(InputError, match="must be identical"):
        layout_from_dict(d)


def test_tee_positions_are_free():
    d = copy.deepcopy(base())
    for b in d["branches"]:
        b["x"] = b["x"] + 100 if b["x"] + 100 < d["main_line"]["length"] else b["x"]
    lay = layout_from_dict(d)
    y = lay.system_y().resolve()
    assert len(y.branch_x) == 13 and y.branch_x.min() == 2200
    x = lay.system_x().resolve()
    assert (x.n_mains_left, x.n_mains_right, x.length) == (6, 9, 10000)


def test_one_sided_layout_puts_the_tee_at_the_end():
    lay = load_layout(LAYOUTS / "M62.yaml")
    x = lay.system_x().resolve()
    assert x.length == 52000 and x.branch_x.tolist() == [52000] and (x.n_mains_left, x.n_mains_right) == (2, 2)
    np.testing.assert_array_equal(x.hanger_x, 52000 - np.array(lay.branch_line("+y").hangers[::-1]))
