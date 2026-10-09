"""Write the plan layout of each 3D model (``inputs/models3d/<M>.json``) as ``inputs/layouts/<M>.yaml``.

Read from the model: the main line (pipe nodes at y = 0, shifted so that it starts at x = 0), its hangers
(fixed nodes at y = 0) and restraints (Pinching4 springs: transverse = along y, longitudinal = along x),
and each branch (tee, side, length to its last pipe node, hangers and restraints, measured from the tee).

    python validation/layouts_from_3d.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from piperom.inputs import dump_yaml  # noqa: E402
from piperom.verification3d import MODELS3D_DIR, list_models3d  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "inputs" / "layouts"
MAT_LONGITUDINAL, MAT_TRANSVERSE = 1, 2     # MinMax-wrapped Pinching4 of the trapezes in the 3D models
N_PIPES = 3                                 # all 3D models: bundles of 3 pipes


def num(v):
    v = round(float(v), 6)
    return int(v) if v.is_integer() else v


def supports(commands) -> dict[tuple[float, float], dict[int, int]]:
    """(x, y) of each hanger -> {direction: material} of its zeroLength spring."""
    coords = {int(a[0]): tuple(a[1:4]) for c, a in commands if c == "node"}
    fixed = {int(a[0]) for c, a in commands if c == "fix"}
    retained = {int(a[2]): int(a[1]) for c, a in commands if c == "rigidLink"}
    out = {}
    for c, a in commands:
        if c == "element" and a[0] == "zeroLength":
            mats = dict(zip(a[a.index("-dir") + 1:], a[a.index("-mat") + 1:a.index("-dir")]))
            for n in (int(a[2]), int(a[3])):
                r = retained.get(n, n)
                if r in fixed:
                    out[coords[r][:2]] = mats
    return out


def layout(name: str) -> dict:
    d = json.loads((MODELS3D_DIR / f"{name}.json").read_text())
    coords = {int(a[0]): tuple(a[1:4]) for c, a in d["commands"] if c == "node"}
    pipe = [coords[n][:2] for n in d["node_groups"]["X"]]
    sup = supports(d["commands"])
    x0 = min(x for x, y in pipe if y == 0)
    L = max(x for x, y in pipe if y == 0) - x0

    def braced(pts, direction, mat):
        return [p for p in pts if sup[p].get(direction) == mat]

    main_h = sorted((x, y) for x, y in sup if y == 0)
    main = {"length": num(L), "n_pipes": N_PIPES,
            "hangers": {"positions": [num(x - x0) for x, _ in main_h]},
            "braces": {"transverse": {"positions": [num(x - x0) for x, _ in braced(main_h, 2, MAT_TRANSVERSE)]},
                       "longitudinal": {"positions": [num(x - x0) for x, _ in braced(main_h, 1, MAT_LONGITUDINAL)]}}}
    branches = []
    for xt in sorted({x for x, y in sup if y != 0}):
        for side, sgn in (("+y", 1), ("-y", -1)):
            hs = sorted(((x, y) for x, y in sup if x == xt and y * sgn > 0), key=lambda p: abs(p[1]))
            if not hs:
                continue
            length = max(abs(y) for x, y in pipe if x == xt and y * sgn > 0)
            branches.append({"x": num(xt - x0), "side": side, "length": num(length),
                             "hangers": {"positions": [num(abs(y)) for _, y in hs]},
                             "braces": {"transverse": {"positions": [num(abs(y)) for _, y in braced(hs, 1, MAT_TRANSVERSE)]},
                                        "longitudinal": {"positions": [num(abs(y)) for _, y in braced(hs, 2, MAT_LONGITUDINAL)]}}})
    return {"name": name,
            "description": f"Plan layout of archetype {name}, read from inputs/models3d/{name}.json "
                           f"(main line shifted by {num(-x0)} mm to start at x = 0)",
            "main_line": main, "branches": branches}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name in list_models3d():
        lay = layout(name)
        (OUT / f"{name}.yaml").write_text("# Units: N, mm, s\n" + dump_yaml(lay), encoding="utf-8")
        print(f"inputs/layouts/{name}.yaml: main {lay['main_line']['length']} mm, {len(lay['branches'])} branches")


if __name__ == "__main__":
    main()
