# The paper's original code

`code implementation for paper/` holds the scripts and results of the paper, as corrected in commit
c13ed81 (2026-10-06) and on 2026-10-09 (see [legacy_issues.md](legacy_issues.md), E5–E6). `piperom` doesn't read from it; only
`validation/` and the tests do, to check that the engine reproduces it ([validation.md](validation.md)).

| Folder | Content |
|---|---|
| `code_proposed_procedure/` | The full procedure of Section 3 on one example system. `equivalent_static.py` performs the equivalent static procedure (with `Pushover2D/Functions.py`) to find the displaced shape. `NLTHA_SDOF.py` builds the equivalent SDOF from it (`Pushover_SDOF/sdof_from_2d.py`) and runs nonlinear time histories under the 150 floor motions. Results in `Results/`. |
| `Section2/` | The OpenSees model used in the demonstration of Section 2. |
| `Pushover2D/` | The equivalent static procedure and the adaptive pushover based on it, for all archetypes of the paper (`Functions.py` + one driver per archetype and direction), with their results. |
| `Pushover_SDOF/` | The equivalent SDOF of each archetype, read from the 2D pushover (`sdof_from_2d.py`): SDOF pushover (`<tag>_SDOF.py`) and time histories under the 44 records × IM1–IM10 (`<tag>_SDOF_NLTHA.py`), with their results. |
| `Results/` | Peak displacements and displaced shapes of the 3D models (`<M>_biron_DispX/Y.txt`, `…DispShapeX/Y.npy`), and the notebooks comparing them with the SDOFs. |
| `3D_models/` | The full 3D models of all archetypes (`<M>_biron.py`), `run_3d.py` to run them in parallel and `results3d.py` to store their results. |
| `trapeze_model_parameters/` | Pinching4 parameters of the transverse and longitudinal trapeze supports (CSV). Copies are the engine's default trapezes (`inputs/trapezes/`). |

The scripts run from any folder; they read the floor motions from `motions/`. Known problems,
inconsistencies and how `piperom` handles each: [legacy_issues.md](legacy_issues.md).

Where `piperom` replaces each part:

| Paper code | `piperom` |
|---|---|
| `Pushover2D/Functions.py`, `equivalent_static.py` | `static_model`, `pushover` |
| `Pushover_SDOF/sdof_from_2d.py` | `sdof` (derived at a chosen Δc) |
| `Pushover_SDOF/<tag>_SDOF.py` (SDOF pushover) | `timehistory.run_sdof_pushover` |
| `Pushover_SDOF/<tag>_SDOF_NLTHA.py`, `code_proposed_procedure/NLTHA_SDOF.py` | `timehistory` |
| `3D_models/*.py` | `inputs/models3d/*.json` + `verification3d` |
| `trapeze_model_parameters/` | `inputs/trapezes/` + `trapeze` |
