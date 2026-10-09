"""Compare the time-history engines with the scripts' stored results and write ``validation/timehistory_report.md``.

* SDOF (archetypes): every record of ``S4_IM`` at IM1-IM10 against
  ``Pushover_SDOF/Results<tag>_SDOF/<tag>_NLTHA_peak_displacements.csv``, (a) with the SDOF the scripts
  analyse (``sdof_from_2d.py``: values read from the pushover results, 3 decimals) and (b) with the engine's
  own SDOF (``derive_sdof`` with the default settings, which round as the scripts).
* SDOF (proposed procedure): the 150 records of ``S4_150`` against
  ``code_proposed_procedure/Results/NLTHA/SDOF_nT4_nL3/Disp_GM<k>_Floor4.txt``, (a) and (b) as above.
* 3D: each 3D model at IM10 under the first ground motion in both orientations, against
  ``Results/<M>_biron_DispX/Y.txt`` (largest peak over the braced nodes, 2 decimals; 100 = collapse).

    python validation/run_timehistory_validation.py [--only sdof|procedure|3d] [--workers N]
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.legacy import ARCHETYPES  # noqa: E402

OUT = Path(__file__).resolve().parent / "timehistory_report.md"
LEVELS = range(1, 11)
SECTIONS = {k: OUT.with_name(f".timehistory_{k}.md") for k in ("sdof", "procedure", "3d")}


def _models(tag: str | None):
    from piperom.inputs import AnalysisSettings, load_system
    from piperom.sdof import derive_sdof
    from piperom.timehistory import SDOFModel
    from validation.legacy import (ARCHETYPE_DIR, EQUIV_STATIC_EXAMPLE, procedure_sdof, procedure_settings,
                                   script_sdof_model)
    if tag is None:
        p = derive_sdof(load_system(EQUIV_STATIC_EXAMPLE), procedure_settings(), 12.0)
        return script_sdof_model(sd=procedure_sdof()), SDOFModel.from_parameters(p)
    p = derive_sdof(load_system(ARCHETYPE_DIR / f"{tag}.yaml"), AnalysisSettings.from_dict())
    return script_sdof_model(tag), SDOFModel.from_parameters(p)


def sdof_task(args) -> tuple:
    tag, set_name, record, level = args
    from piperom.motions import motion_set
    from piperom.timehistory import TimeHistorySettings, run_sdof_time_history
    motion = motion_set(set_name).load(record, level, 4)
    s = TimeHistorySettings.from_dict()
    out = []
    for model in _models(tag):
        r = run_sdof_time_history(model, motion, s)
        out.append((r.peak_u, r.completed, r.collapsed))
    return tag, record, level, out


def _compare(results, ref) -> dict:
    """results: [(key, [(peak, completed, collapsed) script SDOF, ... engine SDOF])]; ref: key -> same."""
    c = {"n": len(results), "same": 0, "flags": 0, "engine_same": 0, "engine_flags": 0}
    for key, (a, b) in results:
        rp, rc, rk = ref[key]
        for prefix, v in (("", a), ("engine_", b)):
            c[prefix + "same"] += f"{v[0]:.4f}" == f"{rp:.4f}"
            c[prefix + "flags"] += (v[1], v[2]) == (rc, rk)
    return c


def sdof_section(workers: int) -> str:
    from piperom.motions import motion_set
    from validation.legacy import script_sdof_peaks
    ms = motion_set("S4_IM")
    tasks = [(tag, "S4_IM", rec, lv) for tag in ARCHETYPES for lv in LEVELS for rec in ms.records]
    by_tag: dict[str, list] = {tag: [] for tag in ARCHETYPES}
    with ProcessPoolExecutor(workers) as pool:
        for tag, rec, lv, out in pool.map(sdof_task, tasks, chunksize=8):
            by_tag[tag].append(((lv, rec), out))
    lines = ["## SDOF time history: archetypes (`Pushover_SDOF/<tag>_SDOF_NLTHA.py`)", "",
             "All 44 records of `S4_IM` × IM1–IM10 (440 runs per archetype), against "
             "`Pushover_SDOF/Results<tag>_SDOF/<tag>_NLTHA_peak_displacements.csv` (peak |u|, 4 decimals; completed "
             "and collapsed flags).", "",
             "- **Scripts' SDOF:** the engine's time history on the SDOF the script analyses (`sdof_from_2d.py`, "
             "which reads Γ, M_eff and φ from the pushover results with 3 decimals). Tests the time-history port.",
             "- **Engine SDOF:** the full engine chain (`derive_sdof` with the default settings, which read the "
             "SDOF as the scripts: pushover step closest to 12 mm, values rounded to 3 decimals).", "",
             "| SDOF | Runs | Scripts' SDOF: identical peaks | same flags | Engine SDOF: identical peaks | "
             "same flags |", "|---|---|---|---|---|---|"]
    for tag in ARCHETYPES:
        c = _compare(by_tag[tag], script_sdof_peaks(tag))
        lines.append(f"| {tag} | {c['n']} | {c['same']} | {c['flags']} | {c['engine_same']} | {c['engine_flags']} |")
        print(f"sdof {tag}: {c}", flush=True)
    return "\n".join(lines)


def procedure_section(workers: int) -> str:
    from piperom.motions import motion_set
    from validation.legacy import procedure_peaks
    ms = motion_set("S4_150")
    peaks = procedure_peaks()
    ref = {(None, rec): (peaks[k], True, False) for k, rec in enumerate(ms.records)}
    with ProcessPoolExecutor(workers) as pool:
        results = [((None, rec), out) for _, rec, _, out in
                   pool.map(sdof_task, [(None, "S4_150", rec, None) for rec in ms.records], chunksize=4)]
    c = _compare(results, {k: v for k, v in ref.items()})
    print(f"procedure: {c}", flush=True)
    return "\n".join([
        "## SDOF time history: proposed procedure (`code_proposed_procedure/NLTHA_SDOF.py`)", "",
        "The SDOF of `equivalent_static.py` (Δc = 12 mm) under the 150 records of `S4_150`, floor 4, against "
        "`Results/NLTHA/SDOF_nT4_nL3/Disp_GM<k>_Floor4.txt` (abs. max of the displacement envelope, 4 decimals). "
        "The stored files have no completed/collapsed flags, so only the peaks are compared.", "",
        "| Runs | Scripts' SDOF: identical peaks | Engine SDOF: identical peaks |", "|---|---|---|",
        f"| {c['n']} | {c['same']} | {c['engine_same']} |"])


def model3d_case(args) -> dict:
    name, rot = args
    from piperom.motions import motion_set
    from piperom.verification3d import Verification3DSettings, braced_peaks, load_model3d, run_3d
    from validation.legacy import script_3d_peak, stored_3d_peak
    ms, m = motion_set("S4_IM"), load_model3d(name)
    a, b = ms.pairs()[0]
    rx, ry = (a, b) if rot == 0 else (b, a)
    t = time.time()
    r = run_3d(m, ms.load(rx, 10), ms.load(ry, 10), Verification3DSettings.from_dict())
    out = {"model": name, "x": rx, "y": ry, "completed": r.completed, "collapsed": r.collapsed,
           "seconds": time.time() - t}
    for d, peak in braced_peaks(m, r).items():
        out[f"new_{d}"] = stored_3d_peak(peak, r.completed, r.collapsed)
        out[f"ref_{d}"] = script_3d_peak(name, d, ms.records.index(rx if d == "x" else ry), 10)
    return out


def model3d_section(workers: int) -> str:
    lines = ["## 3D verification (`piperom.verification3d`)", "",
             "Each model at IM10 under the first ground motion in both orientations, with the default settings, "
             "against `Results/<M>_biron_DispX/Y.txt`: largest peak over the braced nodes as the scripts store it "
             "(2 decimals; 100 = collapse, the analysis stopped when a braced node exceeded 60 mm).", "",
             "| Model | x / y records | Collapsed | x: engine / script (mm) | y: engine / script (mm) | Identical | "
             "Run time (s) |", "|---|---|---|---|---|---|---|"]
    cases = [(m, rot) for m in sorted({t[:3] for t in ARCHETYPES}) for rot in (0, 1)]
    with ProcessPoolExecutor(workers) as pool:
        for c in pool.map(model3d_case, cases):
            print("3d", c, flush=True)
            same = all(f"{c[f'new_{d}']:.2f}" == f"{c[f'ref_{d}']:.2f}" for d in "xy")
            lines.append(f"| {c['model']} | {c['x']} / {c['y']} | {c['collapsed']} | "
                         f"{c['new_x']:.2f} / {c['ref_x']:.2f} | {c['new_y']:.2f} / {c['ref_y']:.2f} | "
                         f"{'yes' if same else 'NO'} | {c['seconds']:.0f} |")
    return "\n".join(lines)


def write_report() -> None:
    parts = ["# Validation of the time-history engines against the scripts", "",
             "Generated by `python validation/run_timehistory_validation.py`. Paths of the scripts are relative to "
             "`code implementation for paper/`; the reference is their state in this commit (c13ed81 and the corrections of 2026-10-09, `docs/legacy_issues.md` E5–E6)."]
    for key, path in SECTIONS.items():
        parts += ["", path.read_text(encoding="utf-8") if path.exists() else f"## {key}: not run yet"]
    OUT.write_text("\n".join(parts) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=list(SECTIONS), help="run one part (the others keep their last result)")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    for key, fn in (("sdof", sdof_section), ("procedure", procedure_section), ("3d", model3d_section)):
        if args.only in (None, key):
            SECTIONS[key].write_text(fn(args.workers), encoding="utf-8")
            write_report()


if __name__ == "__main__":
    main()
