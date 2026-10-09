# `piperom/sdof.py`: equivalent SDOF parameters

## Responsibility

Derives the equivalent SDOF at a chosen Δc: Γ, effective mass, the supports with their shape value φ,
and the scaled Pinching4 spring of each support. It reproduces `Pushover_SDOF/sdof_from_2d.py`
(`load_sdof_params`), which builds the SDOFs of the archetypes and of `code_proposed_procedure/`.

## Interface

**`derive_sdof(system, settings, delta_c=None) -> SDOFParameters`**: runs `pushover.run_step` at
`sdof_delta_c(settings, delta_c)` from a uniform shape.

**`sdof_delta_c(settings, delta_c=None) -> float`**: the requested Δc (default `settings.sdof_delta_c`) or,
with `sdof.on_pushover_grid` (default), the pushover step closest to it. The scripts read the SDOF from
the pushover results at the step closest to 12 mm, i.e. 12.32 mm with the default steps. Each step
starts from a uniform shape, so one step gives the same result as the full pushover.

**`sdof_from_step(rs, step) -> SDOFParameters`**: the same from an existing pushover step, without
re-analysing.

**`SDOFSupport`**:

| Field | Content |
|---|---|
| `kind` | `"transverse"` (braced hanger) or `"longitudinal"` (branch) |
| `dof`, `x`, `branch` | position in the DOF list, coordinate, branch index |
| `phi` | shape value (normalised to the reference DOF); must be > 1e-3 |
| `n_trapezes` | trapezes the spring represents: the number of mains on the hanger's side of `x_center` (transverse), α·`n_braces` (branch) |
| `spring` | the scaled `Pinching4`: deformations / (Γ·φ), forces × `n_trapezes` |

**`SDOFParameters`**:
- fields: `name`, `delta_c`, `gamma`, `effective_mass`, `total_mass`, `mass_ratio`, `u_sdof`,
  `base_shear`, `converged`, `iterations`, `supports` (transverse first, then branches), `d_norm`;
- `n_transverse`, `n_longitudinal_trapezes`, `support_shape` (the `DispShape` vector of the original
  scripts);
- `to_dict()` / `to_yaml()`, `springs_csv()`, `write(out_dir)` (`sdof_parameters.yaml`, `sdof_springs.csv`).

## Design notes

- The spring scales make the SDOF base shear equal to the base shear of the static model, which counts
  the supports of all lumped mains.
- The scripts read Δc, Γ, M_eff and φ from `pushover_results_<tag>.txt`, written with 3 decimals. By
  default the engine rounds them the same way (`sdof.round_decimals`, passed to `sdof_from_step(...,
  decimals=)`); the hysteretic response can be sensitive to it ([legacy_issues.md](../legacy_issues.md), E3).
- The SDOF isn't analysed here; `timehistory.SDOFModel.from_parameters` turns the parameters into
  an analysable model.

## Tests

`tests/test_regression.py`:
- `test_sdof_reproduces_sdof_from_2d`: all 18 archetypes against `load_sdof_params` (Δc, Γ, M_eff, φ and
  scale of every spring, identical);
- `test_sdof_reproduces_the_proposed_procedure_sdof`: against the SDOF `NLTHA_SDOF.py` analysed
  (`sdof_params.json`, identical);
- `test_sdof_reproduces_equivalent_static_procedure`: at full precision against the stored output of
  `code_proposed_procedure/equivalent_static.py` (Γ, M_eff, u_SDOF, shapes and every scaled envelope
  point, relative tolerance 1e-12).
