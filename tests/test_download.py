"""download_floor_motions against a local fake of the Zenodo record (file:// URLs, same JSON layout)."""

import hashlib
import json
import zipfile

import pytest

from piperom.inputs import InputError
from piperom.motions import download_floor_motions


def _fake_record(tmp_path, corrupt=False):
    src = tmp_path / "zenodo"
    src.mkdir(parents=True)
    files = []
    for key, names in (("floor_motions_S4_IM.zip", ["FloorAcc_IM1_120111.txt", "FloorAcc_IM1_120112.txt"]),
                       ("floor_motions_S4_150.zip", ["FloorAcc_1.txt"])):
        zpath = src / key
        with zipfile.ZipFile(zpath, "w") as z:
            for n in names:
                z.writestr("floor_motions/ResultsS4/" + n, "0.001 1 2 3 4\n")
            z.writestr("../evil.txt", "outside")      # must never be extracted
        md5 = hashlib.md5(zpath.read_bytes()).hexdigest()
        files.append({"key": key, "size": zpath.stat().st_size, "checksum": "md5:" + ("0" * 32 if corrupt else md5),
                      "links": {"self": zpath.as_uri()}})
    api = src / "record.json"
    api.write_text(json.dumps({"files": files}))
    return api.as_uri()


def test_download_extracts_and_keeps_existing(tmp_path):
    api = _fake_record(tmp_path)
    dest = tmp_path / "motions"
    existing = dest / "floor_motions" / "ResultsS4" / "FloorAcc_1.txt"
    existing.parent.mkdir(parents=True)
    existing.write_text("original")

    written = download_floor_motions(dest=dest, api_url=api, progress=lambda *_: None)

    assert sorted(p.name for p in written) == ["FloorAcc_IM1_120111.txt", "FloorAcc_IM1_120112.txt"]
    assert existing.read_text() == "original"              # never overwritten
    assert not (tmp_path / "evil.txt").exists()
    assert not list(dest.glob("*.zip"))                     # zips removed by default


def test_download_one_set_and_bad_checksum(tmp_path):
    dest = tmp_path / "motions"
    written = download_floor_motions(["S4_150"], dest=dest, api_url=_fake_record(tmp_path), keep_zip=True,
                                     progress=lambda *_: None)
    assert [p.name for p in written] == ["FloorAcc_1.txt"]
    assert (dest / "floor_motions_S4_150.zip").exists()

    bad = tmp_path / "bad"
    bad.mkdir()
    with pytest.raises(InputError, match="Checksum"):
        download_floor_motions(["S4_IM"], dest=tmp_path / "m2", api_url=_fake_record(bad, corrupt=True),
                               progress=lambda *_: None)
    with pytest.raises(InputError, match="Unknown"):
        download_floor_motions(["S9"], dest=tmp_path / "m3", api_url=_fake_record(tmp_path / "x"))


def test_download_reports_bytes(tmp_path):
    calls = []
    download_floor_motions(["S4_150"], dest=tmp_path / "motions", api_url=_fake_record(tmp_path),
                           progress=lambda *_: None, on_bytes=lambda done, total: calls.append((done, total)))
    assert calls and calls[-1][0] == calls[-1][1]


def test_missing_files(tmp_path, monkeypatch):
    import piperom.motions as mo
    ms = mo.motion_set("S4_IM")
    monkeypatch.setattr(mo, "MOTIONS_DIR", tmp_path)
    present = tmp_path / ms.file_pattern.format(record=ms.records[0], level=1)
    present.parent.mkdir(parents=True)
    present.write_text("0.001 1 2 3 4\n")
    assert len(mo.missing_files(ms, levels=[1])) == len(ms.records) - 1
    assert len(mo.missing_files(ms)) == len(ms.records) * len(ms.levels) - 1
