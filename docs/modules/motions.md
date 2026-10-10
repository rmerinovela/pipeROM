# `piperom/motions.py`: floor-motion sets

## Responsibility

The only access path to floor acceleration histories. It reads the catalogue `motions/motion_sets.yaml`
and loads individual records, converted to mm/s². Format: [input_files.md](../input_files.md#floor-motions).

## Interface

**`MotionSet`**: one catalogue entry.

| Member | Content |
|---|---|
| `name`, `description` | identifier and free text |
| `file_pattern` | relative to `motions/`, with `{record}` and `{level}` |
| `records` | record IDs as strings |
| `record_pairs` | whether consecutive records are the two components of one ground motion |
| `levels` | intensity levels, or `None` |
| `time_column`, `floors` | layout: time column; floor → column |
| `to_mm_s2` | unit factor |
| `path(record, level)` | the file path |
| `pairs()` | `[(x_record, y_record), …]` |
| `load(record, level, floor) -> FloorMotion` | reads one record |

**`FloorMotion`**:
- fields: `set_name`, `record`, `level`, `floor`, `dt` (first value of the time column), `acc` (mm/s²);
- `label`, `duration`.

Functions:
- `load_motion_sets(path=CATALOGUE) -> dict[str, MotionSet]`: an empty dict if there's no catalogue.
- `motion_set(name)`: one set; `InputError` listing the available sets if unknown.
- `select_runs(ms, levels=None, records=None) -> list[(record, level)]`: expands a selection. `None` means
  all; unknown levels or records raise `InputError`.
- `missing_files(ms, levels=None, records=None) -> list[Path]`: files of a selection (default: the whole set)
  that aren't in `motions/`.
- `download_floor_motions(sets=None, dest=MOTIONS_DIR, keep_zip=False, api_url=None, progress=print, on_bytes=None) -> list[Path]`:
  downloads the zips of `sets` (keys of `DATASET_ZIPS`: `S4_IM`, `S4_150`; default both) from the Zenodo
  record `ZENODO_RECORD` (latest version; concept DOI [10.5281/zenodo.23283900](https://doi.org/10.5281/zenodo.23283900)) and extracts them into `dest`. `on_bytes(done, total)` reports the download progress (used by the app's
  progress bar). Returns the files written.

Constants: `ZENODO_DOI`, `ZENODO_RECORD`, `ZENODO_API`, `DATASET_ZIPS` (set → zip name in the record).

## Design notes

- Record IDs are strings; numeric IDs read from text files are normalised (`120111.0` → `"120111"`).
- The time step comes from the file (first time value), as the SDOF scripts (`motion_set(...).load`). The 3D
  scripts compute it as `Timesteps.txt` / 10; for the provided sets the two agree (to the last bit except
  0.0003 s, where they differ by ~5e-20 s).
- Loading a 30,000-row file takes about 0.1–0.2 s; nothing is cached, to keep memory flat in batch runs.
- A missing floor-motion file raises `InputError` with the download command.
- `download_floor_motions` uses only the standard library (`urllib`, `zipfile`, `hashlib`). It reads the
  file list and MD5 checksums from Zenodo's public records API, streams each zip to a `.part` file,
  checks the checksum, and only then extracts. Existing files are never overwritten (the local
  `floor_motions/ResultsS4` may be a link to the original results). Only members under `floor_motions/`
  are extracted, and paths with `..` are rejected.

## Tests

`tests/test_timehistory.py::test_motion_catalogue_and_selection`; `tests/test_download.py` (download
against a local fake of the Zenodo record: extraction, existing files kept, bad checksum, unknown set, byte
progress; `missing_files`); `tests/test_app.py::test_app_offers_floor_motion_download`.
