# `piperom/verification3d.py`: 3D verification models

## Responsibility

Loads the paper's 3D models from `inputs/models3d/*.json` and runs them under a bidirectional floor
motion, with the analysis procedure of the `3D_models/*_biron.py` scripts.

## Interface

**`Model3D`**:
- fields: `name`, `description`, `commands` (`[[opensees_command, [args…]], …]`), `node_groups` (`X`, `Y`,
  `Xt`, `Yt`, …), `rom_systems` (`{"x": system, "y": system}`), `fallback_tests` (per fallback attempt:
  `[test during it, main test restored after it]`, as in the model's script);
- `node_coords()`: node tag → (x, y, z).

Model functions:
- `list_models3d()`, `load_model3d(name)`: `InputError` if missing.
- `model3d_for_system(system_name)`: the 3D model whose `rom_systems` contains the system, or `None`.

**`Verification3DSettings.from_dict(section)`**: gravity steps; damping ratio, `rayleigh_modes` (one mode:
2ξω on the mass term and 2ξ/ω on the stiffness terms; two modes: classical Rayleigh) and factors; solver
`system`; main and fallback convergence tests (`fallback_test: script` uses the model's own types);
`collapse_displacement`, `check_steps` and `record_dt`.

**`build_model3d(model)`**: `op.wipe()`, silences OpenSees, and replays the recorded commands (skipping
`wipe`).

**`run_3d(model, motion_x, motion_y, settings, progress=None) -> Response3D`**:

1. Builds the model.
2. Gravity: load control in `gravity_steps` steps, then `loadConst`.
3. `eigen(4)` for the reported periods; damping from the first mode (default: mass-proportional, 1%).
4. Path series of both components; uniform excitation in x and y.
5. Recorders of the `X` / `Y` node groups every `record_dt` (0.001 s, as the scripts).
6. Newmark transient analysis (`UmfPack`) in chunks of `check_steps` steps; it stops (`collapsed`) when a
   braced node (`Xt` in x, `Yt` in y) exceeds `collapse_displacement`.
7. On failure, the scripts' fallbacks (Newton initial and Broyden dt/2, line search dt/10, Krylov dt/20),
   with each attempt's test types from `fallback_tests`.

**`Response3D`**:
- fields: `model`, `record_x`, `record_y`, `level`, `periods`, `time`, `nodes_x`, `nodes_y`, `ux` and `uy`
  (time × node, relative to the floor), `completed`, `collapsed`, `end_time`, `duration`;
- `peak_x`, `peak_y`: peak |u| per node.

**`braced_peaks(model, response) -> {"x": float, "y": float}`**: the largest peak over the braced nodes
(`Xt` / `Yt`) per direction: the measure stored in `Results/*_DispX/Y.txt`, where the scripts write 100
for a collapsed run and 0 for a run that didn't reach the end (`validation/legacy.stored_3d_peak`).

## Design notes

- **Models as data.** The engine never executes the scripts. `validation/convert_3d_models.py` runs them
  once with a recording stand-in for OpenSees and stores the model-building commands verbatim, with the
  fallback test types parsed from the script, so the conversion can be re-run if a script changes.
- **Fallback tests.** They differ between the scripts ([legacy_issues.md](../legacy_issues.md), E1); the
  engine reproduces each script by default.
- **Cost.** One run takes about 45 s (M01–M03) to several minutes (M29–M63); histories of all
  recorded nodes are kept in memory (tens of MB).

## Tests

`tests/test_timehistory.py::test_3d_verification_reproduces_the_scripts` (M01, IM10, both directions);
`validation/run_timehistory_validation.py` (all models).
