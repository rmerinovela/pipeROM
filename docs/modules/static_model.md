# `piperom/static_model.py`: equivalent static analysis (OpenSees)

## Responsibility

The numerical core: for one assumed shape and one target displacement, it builds the reduced model in
OpenSees, computes the inertial load pattern, solves, and returns the displacements and modal
quantities. It also holds the fixed-point iteration on the shape. Method:
[methodology.md](../methodology.md#1-equivalent-static-procedure-at-a-target-displacement-δc).

It's a restructured port of `build_model` / `iterate_shape_from_static` of the paper's
`Pushover2D/Functions.py`, with:
- the **same OpenSees commands, tags, creation order and floating-point operation order**;
- the regression tests reproducing the paper's results exactly.

Don't "simplify" the arithmetic without re-running `tests/test_regression.py`.

## Interface

**`solve_static_step(rs, d, delta, solver) -> StaticStep`**: one model build and linear solve.

`StaticStep` fields:

| Field | Content |
|---|---|
| `u` | displacements at the DOFs |
| `dof_x`, `dof_mass`, `branch_mass` | DOF positions and masses |
| `loads`, `branch_stay_load` | applied loads |
| `brace_stiffness`, `branch_stiffness` | secant stiffnesses used |
| `gamma`, `effective_mass`, `total_mass`, `mass_ratio` | modal quantities |
| `base_shear`, `f_push` | base shear and load pattern |

`SolverSettings(test, tolerance, max_iterations, branch_split)`: settings of each static step.

**`iterate_shape(rs, d_init, delta, max_iterations, tolerance, solver) -> ShapeResult`**: repeats the solve,
normalising the displacements to max |u| = 1, until the RMS change is below `tolerance`.

`ShapeResult` fields:

| Field | Content |
|---|---|
| `d_star` | converged shape |
| `step` | the last `StaticStep` |
| `iterations`, `converged` | iteration count and status |
| `history` | RMS change per iteration |

Other members:
- `SolverSettings(test, tolerance, max_iterations)`: the OpenSees convergence test of each solve.
- `normalize_shape(phi)`: divides by max |φ|.

## Internal structure of `solve_static_step`

1. Section properties (`_section`) of one main and mass per length; bending stiffness × `n_mains_left` /
   `n_mains_right` on each side of `x_center`.
2. Main-line mesh: nodes at 0, hangers, branch junctions and L. A branch closer than `SNAP_TOL` to a node
   snaps onto it. A massless node is added at `x_center` unless a node is within `SNAP_TOL`.
3. Lumped masses from tributary lengths between mass nodes (hangers, junctions, ends); elastic beam chain,
   each element with the bending stiffness of the side of its midpoint.
4. Hangers:
   - fixed top node, bottom node rigid-linked to the beam;
   - zeroLength spring between them: brace secant stiffness (braced) or `SOFT_STIFFNESS` (unbraced),
     rigid in z and rotations.
5. Branches: a node offset by `BRANCH_OFFSET`, a beam stub (bending stiffness of the mains at the
   junction), the lumped mass and a spring of stiffness α·(`n_braces` / mains at the junction) ×
   longitudinal secant.
6. Load pattern:
   - branch force split into stay and pass parts, in one loop (`_neighbours` finds the tributary
     window). `solver.branch_split` chooses the rule: `consistent` (shares by mass, junction node in the
     main-line share, parts summing to the branch force; the scripts' rule) or `legacy` (the original code);
   - pass part divided left/right;
   - segment-wise redistribution ∝ m·d between branch junctions.
7. Static solve; Γ, M_eff and mass ratio with the shape normalised to the last DOF (the reference branch)
   and the main-line masses × the mains on each DOF's side; base shear = support reactions at the assumed
   shape (braced hangers, branch springs), each × the mains on its side.

The paper's code also ran a diagnostic `eigen(2)` at every solve, whose result nothing used. It's
omitted: results are unchanged (regression tests) and the pushover runs about 3× faster.

## Side effects and constraints

- Uses OpenSees' **global** model: calls `op.wipe()` first and sends OpenSees output to `/dev/null`. Not
  thread-safe; run concurrent analyses in separate processes (`jobs`).
- Tag scheme (`MAIN_NODE`, `HANGER_*`, `BRANCH_*`, …) limits a model to 299 hangers, checked in
  `inputs`.
- Raises `RuntimeError` if the static solve fails (it's linear, so this signals a malformed model).

## Known behaviours kept from the paper's code

- With `branch_split: legacy` only: the stay and pass parts of a branch force don't sum to the branch force
  ([legacy_issues.md](../legacy_issues.md), A2). The default `consistent` split fixes this.
- Branch DOF displacements are read at the main-line junction node, not at the lumped branch node.

## Tests

- `tests/test_regression.py` reproduces all 18 `pushover_results_*.txt` files through `pushover` (default
  settings).
- `tests/test_static_model.py` checks equilibrium (with the consistent split the applied loads equal the
  spring forces of the analysed main), the base shear (support reactions × mains, every archetype) and the
  lumped mains of M29x.
