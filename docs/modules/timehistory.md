# `piperom/timehistory.py`: SDOF nonlinear time history and pushover

## Responsibility

Analyses the equivalent SDOF under a floor acceleration history and reports the response and the
support demands. It's a port of `Pushover_SDOF/<tag>_SDOF_NLTHA.py` (and `code_proposed_procedure/NLTHA_SDOF.py`,
the same procedure): same model, restraint caps, damping, integrator, collapse check and fallback
sequence. It also runs the SDOF pushover of `Pushover_SDOF/<tag>_SDOF.py`.

## Interface

**`SDOFModel`**: what the analysis needs.

| Field | Content |
|---|---|
| `name`, `mass` | identifier and SDOF mass |
| `springs` | `Pinching4` springs acting in parallel |
| `gamma`, `support_phi`, `support_labels`, `support_kinds` | optional: restraint caps, collapse check and support demands need them |

`SDOFModel.from_parameters(p)` builds it from `SDOFParameters`, longitudinal springs first, as the
scripts do. `braced_factor` = Γ·max φ, the largest braced-node displacement per unit SDOF displacement.
Validation builds the scripts' own SDOFs with `validation/legacy.script_sdof_model`.

**`TimeHistorySettings.from_dict(section)`**: the `sdof_time_history` settings merged with the defaults:
damping ratio and Rayleigh factors, convergence test and fallback parameters, `restraint_caps`
(longitudinal, transverse; mm at the support), `collapse_displacement` and `check_steps`.

**`run_sdof_time_history(model, motion, settings) -> SDOFResponse`**:

1. Builds the model:
   - nodes 1 (fixed) and 2 (mass in x and y);
   - one Pinching4 per spring, each wrapped in a MinMax at its cap / (Γ·φ), combined into a `Parallel`
     material;
   - a zeroLength element, rigid in y and rotation.
2. `eigen(1)` → ω; damping coefficients from the Rayleigh factors and ω (default: mass-proportional,
   2ξω with ξ = 1%).
3. Path time series of the floor acceleration (`-prependZero`) and a uniform excitation in x.
4. Newmark (0.5, 0.25) transient analysis (`FullGeneral`, `Transformation`, `EnergyIncr` 1e-4/500) in
   chunks of `check_steps` steps; it stops (`collapsed`) when Γ·max φ·|u| exceeds `collapse_displacement`.
5. On failure, the step-by-step fallbacks, as the scripts (fallback test 1e-3/3000):
   - Newton with initial tangent, dt/2;
   - Broyden, dt/2;
   - Newton with line search, dt/10;
   - Krylov–Newton, dt/20.
6. Displacement and spring force recorded to temporary files and read back.

**`SDOFResponse`**:

| Field | Content |
|---|---|
| `motion`, `record`, `level` | which analysis |
| `period` | SDOF period |
| `time`, `u`, `force` | histories; `u` relative to the floor |
| `completed`, `collapsed`, `end_time`, `duration` | analysis status |
| `peak_u` | stored at creation, so it survives when histories are dropped to save memory |

**`run_sdof_pushover(model, targets=None, step=0.05) -> SDOFPushover`**: displacement-controlled pushover of
the SDOF (springs without caps) to 50 mm, or through a list of targets (`CYCLIC_PROTOCOL`, the scripts'
cyclic option), with the scripts' fallbacks (sub-steps of step/2, step/2, step/10). Returns `u`, `force`
(= − support reaction) and `completed`.

Other functions:
- `sdof_period(model) -> float`: initial period.
- `support_demands(model, peak_u) -> list[dict]`: per support, peak displacement Γ·φ·u and its ratios to
  the trapeze envelope deformations ePd1–ePd4.

## Side effects and constraints

- Uses OpenSees' global model (`op.wipe()` at start and end) and a temporary directory for recorders.
  Run concurrent analyses in separate processes.
- About 0.3 s per 30 s floor record (30,000 steps).

## Tests

`tests/test_timehistory.py`: M01x, M03x (collapses) and M62y against the scripts' stored peaks and flags;
the proposed procedure against `NLTHA_SDOF.py`; the SDOF pushover of M01x, M29x and M62y against
`DispC.out` / `VbaseC.out`. `validation/run_timehistory_validation.py` covers every archetype, record and
IM1–IM10, and the 150 records of the proposed procedure.
