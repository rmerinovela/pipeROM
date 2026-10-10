"""Smoke test of the Streamlit app (headless)."""

from pathlib import Path

import pytest

st_testing = pytest.importorskip("streamlit.testing.v1")
APP = str(Path(__file__).resolve().parents[1] / "app" / "streamlit_app.py")
HAVE_MOTIONS = (Path(__file__).resolve().parents[1] / "motions" / "floor_motions" / "ResultsS4").exists()


def button(at, label):
    return next(b for b in at.button if b.label == label)


def test_app_runs_pushover_and_sdof_in_both_directions():
    at = st_testing.AppTest.from_file(APP, default_timeout=240)
    at.run()
    assert not at.exception, at.exception
    assert not at.error, [e.value for e in at.error]

    # switch hanger definitions back and forth: inputs must survive
    at.radio(key="m_hmode").set_value("Explicit positions").run()
    at.radio(key="m_hmode").set_value("Regular grid").run()
    assert not at.exception, at.exception
    assert at.session_state.m_hfirst == 1000

    at.number_input(key="s_n").set_value(4).run()
    button(at, "Run pseudo-pushover").click().run()
    assert not at.exception, at.exception
    assert {d: len(r.steps) for d, r in at.session_state.po_result.items()} == {"x": 4, "y": 4}

    button(at, "Derive SDOF parameters").click().run()
    assert not at.exception, at.exception
    sd = at.session_state.sdof_result
    # M02 layout: y = main line (4 transverse restraints), x = branch line (2 transverse restraints)
    assert (sd["y"].n_transverse, sd["x"].n_transverse) == (4, 2)
    assert sd["y"].delta_c == round(0.1 + 49.9 / 3, 3)    # closest pushover step, rounded as the scripts


def test_app_adds_a_branch_side():
    at = st_testing.AppTest.from_file(APP, default_timeout=120)
    at.run()
    at.text_input(key="n_tees").set_value("18000").run()
    assert not at.exception, at.exception
    # the new -y branch has no restraints yet: a clear message, not a crash
    assert any("no restraint along y" in e.value for e in at.error)
    at.text_input(key="n_tees").set_value("").run()
    assert not at.error, [e.value for e in at.error]


@pytest.mark.skipif(not HAVE_MOTIONS, reason="floor motions not available")
def test_app_runs_time_history_in_both_directions():
    at = st_testing.AppTest.from_file(APP, default_timeout=600)
    at.run()
    button(at, "Derive SDOF parameters").click().run()
    at.toggle(key="th_all").set_value(False).run()
    button(at, "Run SDOF time histories").click().run()
    assert not at.exception, at.exception
    for d in ("x", "y"):
        model, resp = at.session_state.th_result[d]
        assert len(resp) == 2 and all(r.completed or r.collapsed for r in resp)


def test_app_reports_invalid_layout():
    at = st_testing.AppTest.from_file(APP, default_timeout=60)
    at.run()
    at.radio(key="m_hmode").set_value("Explicit positions").run()
    at.text_area(key="m_hpos").set_value("1000, 500").run()
    assert not at.exception, at.exception
    assert any("increasing" in e.value for e in at.error)
    at.text_input(key="p_tees").set_value("40000").run()
    assert any("tee must lie on the main line" in e.value or "increasing" in e.value for e in at.error)


@pytest.mark.skipif(not HAVE_MOTIONS, reason="floor motions not available")
def test_app_offers_floor_motion_download():
    at = st_testing.AppTest.from_file(APP, default_timeout=240)
    at.run()
    button(at, "Derive SDOF parameters").click().run()
    assert not at.exception, at.exception
    dl = [b for b in at.button if b.label.startswith("Download S4_IM from Zenodo")]
    assert {b.key for b in dl} == {"dl_motions_th", "dl_motions_3d"}
    assert all(b.disabled for b in dl)          # every file present locally: nothing to download
    assert any("available locally" in c.value for c in at.caption)
