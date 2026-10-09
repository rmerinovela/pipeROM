# Known issues in the original scripts

The scripts are in `code implementation for paper/`; paths below are relative to that folder. This file lists
what was found while reimplementing them in `piperom`, and where each issue stands:

- **Scripts (c13ed81):** whether the scripts were corrected in commit c13ed81 (2026-10-06), which is the
  reference the engine now reproduces ([validation.md](validation.md)).
- **piperom:** how the engine handles it.

Line numbers refer to the scripts before c13ed81.

## A. Errors that change results

| # | Where | Issue | Scripts (c13ed81) | piperom |
|---|---|---|---|---|
| A1 | `Pushover2D/Functions.py:43`, `code_proposed_procedure/equivalent_static.py:87` | The static procedure's longitudinal trilinear used **11500 N at 24 mm**; the Pinching4 parameters (CSV, SDOF scripts, 3D models) have **10000 N at 24 mm** (11500 N is at 61 mm). | **Fixed**: 10000 N; `equivalent_static.py` now uses `Functions.py`. | The trilinear is derived from the trapeze files (10000 N). |
| A2 | `Functions.py:630-830`, `equivalent_static.py:670-870` | Branch force split: the "stay" share included the main-line junction mass in its denominator, the "pass" share didn't (line 788 commented out), so stay + pass ≠ branch force. | **Fixed** (line uncommented). | Default `equivalent_static.branch_split: consistent` (the fix). `legacy` keeps the original split. |
| A3 | `NLTHA_SDOF.py:193` | Longitudinal strength per branch = `nL / n_ortho` instead of each branch's own `n_braces`. | **Fixed**: `NLTHA_SDOF.py` builds the SDOF with `Pushover_SDOF/sdof_from_2d.py`, per branch. | Each branch's α·`n_braces`. |
| A4 | `Pushover_SDOF/*.py` | Hand-typed SDOF constants not consistently derived from the 2D pushover (Δc not recorded and varying, mass 1–12% low, shapes and nL not matching the 2D model, M62y/M63y ignoring a branch). | **Fixed**: every SDOF is read from the pushover results (`sdof_from_2d.py`, step closest to Δc = 12 mm) and its layout file. | `piperom sdof` derives the SDOF from the inputs; by default at the pushover step closest to Δc (`sdof.on_pushover_grid`), as the scripts. |
| A5 | `Pushover_SDOF/*.py:238,246,254` | Fallbacks called `op.analyze(1, dsteps/2)` in a static analysis (smaller step ignored) and left the test as `RelativeEnergyIncr` afterwards. | **Fixed**: sub-steps through `DisplacementControl`, main test restored. | `run_sdof_pushover` (`piperom sdof --pushover`), as the corrected scripts. |
| A6 | `3D_models/M01-M03_biron.py` | The main test `RelativeEnergyIncr` (1e-4) never converges from rest; the run stopped at t ≈ 0.002 s. | **Fixed**: main test `EnergyIncr` in all models. The fallback tests still differ between models (E1). | Main test `EnergyIncr`; fallback tests of each model's script (E1). |
| A7 | `Pushover_SDOF/M61*`, `M62*` vs `Results/M61_*`, `M62_*` | The paper's SDOF peaks of M61x/y and M62x/y weren't reproduced by any available SDOF definition. | **Resolved**: all SDOFs and their time histories were regenerated (`Results<tag>_SDOF/*_NLTHA_peak_displacements.csv`). | Reproduced exactly ([validation.md](validation.md)). |

## B. Crashes and dead code paths

| # | Where | Issue | Scripts (c13ed81) | piperom |
|---|---|---|---|---|
| B1 | `NLTHA_SDOF.py:72` | `nT`, `nL` used before definition: `NameError`. | **Fixed** (script rewritten). | n/a |
| B2 | `Functions.py:702-714` | Segment variables defined only inside the branch loop: a system **without branches** crashed. | **Partly**: the variables are now defined outside the loop; `d_all / d_all[-1]` still assumes the last DOF is a branch. | Input error: at least one branch is required. |
| B3 | `Functions.py:448` (now 496) | A branch with zero length, pipes or braces (or `alpha <= 0`) is skipped, leaving `None` in the node list. | Unchanged. | Input validation rejects these values. |
| B4 | `equivalent_static.py:9-44` | `compute_stiff_mask` built its own hanger list without the 1 m end clearance. | **Fixed**: same hanger rule as `build_model`. | Uses the model's hanger list. |
| B5 | `Functions.py` (tags) | Node/element tags overlap above ~300 hangers. | Unchanged. | Limit checked (299 hangers). |
| B6 | `Functions.py` (`iterate_shape_from_static`) | `build_model` is called once before the loop and its result discarded. | Unchanged (no effect on results). | Removed. |
| B7 | `Pushover2D/CS_*.py` | `d_init_current` is updated but `d0` is always passed: every Δc step starts from a uniform shape. | Unchanged. | `warm_start` setting, default `false` as in the scripts. |

## C. Input/output and portability

| # | Where | Issue | Scripts (c13ed81) |
|---|---|---|---|
| C1 | `NLTHA_SDOF.py`, `3D_models/*`, `Results/*.ipynb` | Hard-coded `C:/Users/rmeri/...` paths; inputs read relative to the working directory. | **Fixed** in the scripts (paths relative to each script; motions from `motions/`). One notebook still saves figures to a local folder (`Results/Summary_all_archetypes.ipynb`). |
| C2 | `Pushover2D/CS_*.py` | Output column order differed between drivers. | **Fixed**: each driver writes `pushover_layout_<tag>.json` with the column order. |
| C3 | `equivalent_static.py`, `NLTHA_SDOF.py` | Hand-off through `Results/EquivStatic` relative to the working directory. | **Fixed**: next to the scripts (`code_proposed_procedure/Results/`). |
| C4 | `Functions.py` | Debug prints on every iteration; `equivalent_static.py` duplicated `Functions.py`. | **Fixed** (the iteration log remains). |
| C5 | `Section2/ContinuousSystem_3DOF.py`, `Pushover_SDOF/*_SDOF.py` | `vfo` 0.0.19 fails with numpy ≥ 2 when plotting. | `Pushover_SDOF/*_SDOF.py` still import `vfo` (unused). |
| C6 | `Pushover_SDOF/Results*/*.out` | Windows line endings. | Regenerated. |
| C7 | `Results/*_DispX.txt`, `*_DispY.txt` | Undocumented. | Documented in `3D_models/results3d.py`: largest peak over the braced nodes, row = record applied in that direction, columns IM1–IM10, 100 = collapse, 0 = not completed. |
| C8 | `NLTHA_SDOF.py` | Read 150 records from a folder outside the repository. | **Fixed**: `motion_set('S4_150')`. That `motions/floor_motions/ResultsS4/FloorAcc_1`–`150.txt` is the set of the paper is assumed, not verified. |

## D. Worth checking (possibly intentional)

| # | Where | Observation | Scripts (c13ed81) |
|---|---|---|---|
| D1 | `3D_models/*.py:316-317` | `MinMax` caps looked swapped (longitudinal ±35 mm, transverse ±60 mm). | **Fixed**: longitudinal 60, transverse 35; the SDOF springs have the same caps. |
| D2 | `Results/M01_biron.ipynb` cell 22 | `np.flip(POx[:,6:13])` flips both axes. | Not re-checked. |
| D3 | `Results/M01_biron.ipynb` cells 9 and 11 | Different Γ·φ for x in two plots. | Not re-checked. |
| D4 | `NLTHA_SDOF.py` vs `3D_models` | Different damping models in the SDOF and the 3D benchmark. | **Fixed**: both 1% mass-proportional at the first (elastic) mode. |
| D5 | `NLTHA_SDOF.py:285-290` | Upward "gravity" load on a rigid DOF; "Floor Motion Done" printed even on failure. | **Fixed** (no gravity step in the SDOF scripts). |
| D6 | `Functions.py:93` | `trilinear_keff` returns the secant stiffness; comments said "tangent". | **Fixed** (comments). |

## E. Found when porting the c13ed81 scripts (2026-10-09)

| # | Where | Observation | piperom |
|---|---|---|---|
| E1 | `3D_models/*_biron.py`, fallback loop | The convergence-test types of the four fallback attempts differ between models: M01–M03 use `RelativeEnergyIncr` throughout; M29–M31 use `EnergyIncr` during the attempt and then set the main test to `RelativeEnergyIncr`; M61–M63 use `RelativeEnergyIncr` during the first three and `EnergyIncr` after them, the reverse for `KrylovNewton`. After a fallback, the main test can stay `RelativeEnergyIncr` for the rest of the run. They only act when a step fails. | Stored per model (`fallback_tests` in `inputs/models3d/*.json`, from `validation/convert_3d_models.py`) and used by default (`verification_3d.fallback_test: script`); a single test type can be set instead. |
| E2 | `3D_models/*_biron.py` | Displacements are recorded every 0.001 s. For the records with a time step of 0.0003 or 0.0005 s the peak is taken from these samples. | `verification_3d.record_dt: 0.001` (null = every step). |
| E3 | `Pushover_SDOF/sdof_from_2d.py` | The SDOF is read from `pushover_results_<tag>.txt`, written with 3 decimals (Γ, M_eff, φ). The rounding is small (e.g. Γ 0.90185 → 0.902), but the hysteretic response can switch path within such an interval: for M01x, IM9, record 120522 the peak goes from 23.86 mm (full precision) to 20.94 mm (rounded Γ), while ±0.01% around either value changes it by < 0.01%. | Reproduced: by default the engine rounds the SDOF the same way (`sdof.round_decimals: 3`; null = full precision). With full precision, most peaks differ slightly and a few records by up to ~14%, with the same collapse flags. |
| E4 | `Pushover_SDOF/ResultsM61y_SDOF/DispC.out`, `VbaseC.out`, `pushover_comparison_M61y.png` | Written on 2026-10-02, before `Pushover2D/pushover_results_M61y.txt` was regenerated (2026-10-06): the stored SDOF pushover of M61y doesn't match its current SDOF (up to 3.4% in force). Re-running `M61y_SDOF.py` gives the engine's curve exactly. | **Fixed 2026-10-10**: `M61y_SDOF.py` rerun; its files now match the engine's curve to the 6 digits stored. |
| E5 | `Pushover2D/CS_Lumped_Iter_M29*-M31*_cont_PO.py` (c13ed81) | The 2D models of M29–M31 didn't match their 3D models: branches 5.5 m instead of 5 m; the line analysed in x 12 m (hangers at 1.5 … 10.5 m) instead of the 10 m between the branch ends; the main line 31.5 m on a regular hanger grid instead of 32 m with the 3D hangers (0, 2, 5 … 20, 23.125, 26, 29 m); M30y/M31y in 3D coordinates shifted by 1.05 m; the two branches sharing a tee at 7.05 and 11.55 m as one branch of 3 pipes (M29y, half their mass) or as two DOFs at the same point (M30y, M31y). | **Fixed 2026-10-09**: the drivers follow the 3D models (`Functions.build_model` accepts explicit hanger positions, `x_hangers_user`); the shared tees are one branch of 6 pipes and 2 braces; pushovers, SDOFs, SDOF pushovers and SDOF time histories of M29–M31 regenerated. The archetype files follow (`hangers: positions`). |
| E6 | `Pushover2D/CS_Lumped_Iter_M0*x`, `M6*x`, `M61y`; `3D_models/M01_biron.py`, `M61_biron.py` | 2D models inconsistent with their 3D models (found 2026-10-09): in the x models of M01–M03 and M61–M63 the hangers and transverse braces along the branch were measured from its free end instead of the tee (mirrored); M62x's branch braces also differed; M61x used 8 of the main line's 9 longitudinal braces; M61y's main-line braces differed in 6 of 10 positions. | **Fixed 2026-10-09.** M02x, M03x, M61x, M62x, M63x: the 2D drivers follow the 3D models (2D pushovers, SDOFs and SDOF time histories regenerated). M01x, M61y: the 3D models follow the roughly evenly spaced 2D models instead, as the design procedure will space restraints evenly (M01: branch hangers at 2, 5 … 17 m from the tee; M61: main-line y restraints at 4, 10, 19, 28, 37, 46, 52, 61, 70, 79 m); their 3D NLTHA rerun (all 220 runs each; reproduced by the engine). All 18 directions now reduce from the plan layouts (`inputs/layouts/`) to the paper's 2D models. Final RE medians: M01x 1.16 (was 1.19), M61y 1.19 (was 1.14), M01y 0.99, M61x 0.95; overall, median RE over all intensities (collapse included) within 10% for 15 of 18 directions and within 20% for 18 of 18 (132 and 166 of 180 points); excluding the intensities at which ≥ 50% of the 3D records collapsed, 14 and 18 of 18 (120 and 154 of 160 points). At 6 of the 20 collapse intensities the SDOF doesn't collapse (RE 0.30–0.48: M03, M62, M63). Why consistent layouts can fit worse: the converged equivalent static shape (load ∝ m·d, i.e. the first mode at Δc) amplifies long, lightly braced spans and unrestrained overhangs (max φ 1.20 after one iteration, 1.62–1.66 converged, for the inconsistent alternatives tested; 3D static 1.00–1.08, 3D NLTHA 1.0–1.2); Γ partly absorbs this amplitude error in the demand Γφu (shape form and amplitude per archetype: `Results/Summary_all_archetypes.ipynb`). |
