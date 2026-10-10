"""Floor-motion sets (read only from the ``motions/`` folder, described in ``motions/motion_sets.yaml``)."""

from __future__ import annotations

import hashlib
import json
import shutil
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from . import REPO_DIR
from .inputs import InputError

MOTIONS_DIR = REPO_DIR / "motions"
CATALOGUE = MOTIONS_DIR / "motion_sets.yaml"

# Zenodo dataset with the floor-motion files (too large for git). ZENODO_DOI is the concept DOI (all
# versions); the download uses the latest version of the record.
ZENODO_DOI = "10.5281/zenodo.23283900"
ZENODO_RECORD = "23283901"     # any version of the record; its latest version is downloaded
ZENODO_API = "https://zenodo.org/api/records/"
DATASET_ZIPS = {"S4_IM": "floor_motions_S4_IM.zip", "S4_150": "floor_motions_S4_150.zip"}


@dataclass
class FloorMotion:
    set_name: str
    record: str
    level: int | None
    floor: int
    dt: float
    acc: np.ndarray            # mm/s^2, first value at time dt (OpenSees Path series with -prependZero)

    @property
    def label(self) -> str:
        lv = "" if self.level is None else f"IM{self.level} "
        return f"{lv}{self.record} (floor {self.floor})"

    @property
    def duration(self) -> float:
        return self.dt * len(self.acc)


@dataclass
class MotionSet:
    name: str
    description: str
    file_pattern: str
    records: list[str]
    record_pairs: bool
    levels: list[int] | None
    time_column: int
    floors: dict[int, int]
    to_mm_s2: float

    def path(self, record: str, level: int | None = None) -> Path:
        if self.levels is not None and level is None:
            raise InputError(f"Motion set '{self.name}' needs an intensity level")
        return MOTIONS_DIR / self.file_pattern.format(record=record, level=level)

    def pairs(self) -> list[tuple[str, str]]:
        if not self.record_pairs:
            raise InputError(f"Motion set '{self.name}' has no record pairs")
        return [(self.records[i], self.records[i + 1]) for i in range(0, len(self.records) - 1, 2)]

    def load(self, record: str, level: int | None = None, floor: int = 4) -> FloorMotion:
        if record not in self.records:
            raise InputError(f"Record '{record}' is not in motion set '{self.name}'")
        if self.levels is not None and level not in self.levels:
            raise InputError(f"Level {level} is not available in motion set '{self.name}'")
        if floor not in self.floors:
            raise InputError(f"Floor {floor} is not available in motion set '{self.name}'")
        path = self.path(record, level)
        if not path.exists():
            raise InputError(f"Missing floor-motion file: {path.relative_to(REPO_DIR)}. "
                             f"Download the floor motions with: python -m piperom download-motions")
        data = np.loadtxt(path)
        return FloorMotion(self.name, record, level, floor, dt=float(data[0, self.time_column]),
                           acc=self.to_mm_s2 * data[:, self.floors[floor]])


def _records(spec, name: str) -> list[str]:
    if isinstance(spec, list):
        return [str(r) for r in spec]
    if isinstance(spec, dict) and "file" in spec:
        text = (MOTIONS_DIR / spec["file"]).read_text().split()
        return [str(int(float(v))) if float(v).is_integer() else v for v in text]
    if isinstance(spec, dict) and "range" in spec:
        a, b = spec["range"]
        return [str(i) for i in range(int(a), int(b) + 1)]
    raise InputError(f"Motion set '{name}': 'records' must be a list, {{file: ...}} or {{range: [a, b]}}")


def load_motion_sets(path: Path = CATALOGUE) -> dict[str, MotionSet]:
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text()) or {}
    out = {}
    for name, d in data.items():
        try:
            out[name] = MotionSet(
                name=name, description=" ".join(str(d.get("description", "")).split()),
                file_pattern=d["file_pattern"], records=_records(d["records"], name),
                record_pairs=bool(d.get("record_pairs", False)),
                levels=None if d.get("levels") is None else [int(v) for v in d["levels"]],
                time_column=int(d.get("time_column", 0)),
                floors={int(k): int(v) for k, v in d["floors"].items()},
                to_mm_s2=float(d["to_mm_s2"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise InputError(f"Motion set '{name}' in {path.name} is invalid: {exc}") from None
    return out


def motion_set(name: str) -> MotionSet:
    sets = load_motion_sets()
    if name not in sets:
        raise InputError(f"Unknown motion set '{name}'. Available: {', '.join(sets) or 'none'}")
    return sets[name]


def select_runs(ms: MotionSet, levels: list[int] | None = None,
                records: list[str] | None = None) -> list[tuple[str, int | None]]:
    """(record, level) pairs to analyse; ``None`` means every record / every level of the set."""
    if ms.levels is None:
        lv: list[int | None] = [None]
    else:
        lv = list(levels) if levels else list(ms.levels)
        bad = [v for v in lv if v not in ms.levels]
        if bad:
            raise InputError(f"Levels {bad} are not available in motion set '{ms.name}'")
    recs = [str(r) for r in records] if records else list(ms.records)
    bad = [r for r in recs if r not in ms.records]
    if bad:
        raise InputError(f"Records {bad} are not in motion set '{ms.name}'")
    return [(r, l) for l in lv for r in recs]


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def missing_files(ms: MotionSet, levels: list[int] | None = None, records: list[str] | None = None) -> list[Path]:
    """Floor-motion files of a selection (default: the whole set) that are not in ``motions/``."""
    return [p for p in (ms.path(r, lv) for r, lv in select_runs(ms, levels, records)) if not p.exists()]


def download_floor_motions(sets: list[str] | None = None, dest: Path = MOTIONS_DIR, keep_zip: bool = False,
                           api_url: str | None = None, progress=print, on_bytes=None) -> list[Path]:
    """Downloads the floor-motion files of ``sets`` (default: all) from the Zenodo dataset (``ZENODO_DOI``)
    (latest version) and extracts them into ``dest`` (``motions/``). Each zip is checked against the MD5 checksum published
    by Zenodo. Files that already exist are kept, never overwritten. ``on_bytes(done, total)``, if given, is
    called while each zip downloads. Returns the files written."""
    names = list(DATASET_ZIPS) if not sets else list(sets)
    bad = [s for s in names if s not in DATASET_ZIPS]
    if bad:
        raise InputError(f"Unknown floor-motion sets {bad}. Available: {', '.join(DATASET_ZIPS)}")

    url = api_url or f"{ZENODO_API}{ZENODO_RECORD}/versions/latest"
    try:
        with urllib.request.urlopen(url) as r:
            record = json.load(r)
    except OSError as exc:
        raise InputError(f"Could not read the Zenodo record {url}: {exc}") from None
    files = {f["key"]: f for f in record.get("files", [])}

    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    written = []
    for s in names:
        key = DATASET_ZIPS[s]
        if key not in files:
            raise InputError(f"{key} is not in the Zenodo record {url}")
        info = files[key]
        zpath = dest / key
        if not (zpath.exists() and _md5(zpath) == info["checksum"].split(":")[-1]):
            progress(f"downloading {key} ({info['size'] / 1e6:.0f} MB) from doi:{ZENODO_DOI} ...")
            part = zpath.with_name(key + ".part")
            with urllib.request.urlopen(info["links"]["self"]) as r, open(part, "wb") as f:
                done = 0
                while chunk := r.read(1 << 20):
                    f.write(chunk)
                    done += len(chunk)
                    if on_bytes is not None:
                        on_bytes(done, info["size"])
            if _md5(part) != info["checksum"].split(":")[-1]:
                part.unlink()
                raise InputError(f"Checksum mismatch for {key}: the download is corrupt, try again")
            part.replace(zpath)

        new = skipped = 0
        with zipfile.ZipFile(zpath) as z:
            for m in z.infolist():
                target = dest / m.filename
                if m.is_dir() or not m.filename.startswith("floor_motions/") or ".." in Path(m.filename).parts:
                    continue
                if target.exists():
                    skipped += 1
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(m) as src, open(target, "wb") as out:
                    shutil.copyfileobj(src, out, 1 << 20)
                written.append(target)
                new += 1
        progress(f"{s}: {new} files extracted, {skipped} already present (kept)")
        if not keep_zip:
            zpath.unlink()
    return written
