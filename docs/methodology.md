# Methodology: what the engine computes

This file describes the analyses `piperom` can run, the assumptions they keep from the paper, and their
limits. Units are N, mm, s; masses are in tonnes.

The engine implements the reduced-order model of Merino, Gentile and Galasso, "A reduced-order model
for performance-based seismic design of suspended piping systems in buildings". It follows the paper's
scripts (`code implementation for paper/`, as of commit c13ed81 with the corrections of 2026-10-09) step by step and reproduces their stored
results (see the end of this file); only inputs and outputs differ.

```
system file ──► 1. equivalent static procedure ──► 2. adaptive pseudo-pushover
(one direction)        (at one Δc)                        (Δc = Δc1 … Δcn)
                           │
                           ▼
                  3. equivalent SDOF at a chosen Δc ──► 4. SDOF time history under floor motions
                                                                       │
            5. 3D verification (paper archetypes only) ◄── compared ───┘
```

## The idealised system

A system file describes **one loading direction**:

- **Main line.** A straight pipe bundle (`n_pipes` identical pipes) along x, loaded transversely (y). It's
  modelled as an elastic beam with lumped masses at its nodes:
  - mass per unit length = `mass_factor` × (pipe + contained fluid);
  - defaults: factor 1.35, water = steel density / 7.8;
  - nodal masses from tributary lengths between mass nodes (hangers, branch junctions, ends).
- **Lumped mains.** `n_mains` identical parallel mains can be lumped into the analysed one, with a
  different number on each side of `x_center` (`n_mains_left`, `n_mains_right`; default length / 2). The
  analysed main keeps the mass of one main; the others:
  - add their bending stiffness (I × n_mains on each side; a massless node at `x_center`);
  - share each branch's longitudinal braces (a branch spring = `n_braces / n_mains` trapezes);
  - add their mass to Γ, M_eff and the total mass, and their supports to the base shear.

  A point at `x_center` belongs to both sides (the average of the two numbers).
- **Gravity hangers.** At given positions along the main line:
  - no lateral stiffness;
  - rigid vertically and in rotation.
- **Transverse braces (trapezes).** A subset of the hangers. Each one is a spring in the transverse
  direction (and along the pipe axis, as in the paper's code).
- **Branches.** Straight pipelines orthogonal to the main line, framing into it at a junction and loaded
  along their axis. Each branch is lumped into **one DOF**:
  - its whole mass;
  - `n_braces` longitudinal trapezes in parallel.

  A branch can't coincide with a braced hanger. At least one branch is required: the procedure is
  defined for a main line plus branches.
- **Trapeze behaviour.** One Pinching4 definition per type (transverse, longitudinal), shared by all
  braces of that type. The static procedure uses a trilinear secant idealisation of its positive envelope:
  - origin → envelope point 2 → envelope point 3;
  - then a slope of 1% of the first branch.

**DOFs:** all hangers (by position), then all branches (in input order). The **last branch** is the
reference DOF used to normalise the shape. Γ depends on this choice; the SDOF definition (Γ·φ and the
effective mass) does not.

## 1. Equivalent static procedure (at a target displacement Δc)

For an assumed normalised shape `d` (max |d| = 1):

1. **Stiffness.** Every brace and branch spring gets the secant stiffness of its trilinear backbone at
   displacement Δc·dᵢ.
2. **Loads.** The spring forces implied by the shape become an inertial load pattern:
   - Each branch force splits into a part that stays on the branch and a part passed to the main line,
     by tributary mass: the branch mass against the main-line mass around the junction, junction node
     included. The two parts add up to the branch force. The original code left the node mass out of the
     passed part, so the parts didn't add up; `branch_split: legacy` reproduces that
     ([legacy_issues.md](legacy_issues.md), A2).
   - The passed part divides left/right of the junction.
   - Along the main line, between consecutive branch junctions, the total shear is redistributed in
     proportion to mᵢ·dᵢ.
3. **Solve.** A linear static analysis (OpenSees) gives the displaced shape; normalising it gives a new `d`.
4. **Iterate.** Repeat until the RMS change of `d` is below the tolerance (fixed-point iteration).

**Outputs:**
- the converged shape;
- the participation factor Γ = Σmd / Σmd² and the effective mass M_eff = Γ·Σmd (shape normalised to
  the reference DOF; m includes all lumped mains);
- the effective-mass ratio and the base shear V_b = the support reactions at the assumed shape (braced
  hangers and branch springs), each times the number of mains on its side;
- the equivalent SDOF displacement u_SDOF = Δc / (Γ · max d).

## 2. Adaptive pseudo-pushover

The procedure is repeated for a list of target displacements, Δc = 0.1 … 50 mm by default. Each step
gives one point (u_SDOF, V_b) of the capacity curve, with its own shape: the shape adapts as the braces
soften. By default each step starts from a uniform shape, as in the paper's code (`warm_start: false`).

## 3. Equivalent SDOF

At a chosen Δc (default 12 mm), from the converged shape. The scripts read the SDOF from the pushover
results (`Pushover_SDOF/sdof_from_2d.py`), so by default the engine does the same:
- Δc is the pushover step closest to the requested value (`sdof.on_pushover_grid`): 12.32 mm with the
  default steps;
- Δc, Γ, M_eff and the shape are rounded to 3 decimals, as written in those results
  (`sdof.round_decimals`; null keeps full precision).

The rounding matters for some records: the hysteretic response can switch path within it (e.g. Γ 0.90185
vs 0.902 changes one M01x peak by 12%; [legacy_issues.md](legacy_issues.md), E3).

| Quantity | Value |
|---|---|
| Mass | M_eff |
| Springs (in parallel) | one per support: each braced hanger (the number of mains on its side of `x_center`) and each branch (α·`n_braces` trapezes) |
| Spring definition | the support's Pinching4, envelope **deformations divided by Γ·φᵢ**, **forces multiplied by the number of trapezes** |

φᵢ is the shape at the support, normalised to the reference DOF; it must be positive. The SDOF base shear
equals the base shear of the static model, so the SDOF backbone reproduces the capacity curve. A
displacement-controlled SDOF pushover (to 50 mm, or cyclic) checks this, as `Pushover_SDOF/<tag>_SDOF.py`.

## 4. SDOF time history

The SDOF is excited by a floor acceleration history (`motions/`), as in `Pushover_SDOF/<tag>_SDOF_NLTHA.py`:

- a single zeroLength element with the parallel springs;
- each spring wrapped in a MinMax at the 3D model's restraint cap (longitudinal 60 mm, transverse 35 mm
  at the support, i.e. cap / (Γ·φᵢ) in SDOF displacement);
- mass-proportional damping, 1% of critical at the initial (elastic) SDOF frequency (the trapezes
  dissipate through their hysteresis);
- Newmark average-acceleration integration, in chunks of 100 steps: the analysis stops (collapse) when
  the largest braced-node displacement Γ·max φ·u exceeds 60 mm;
- the scripts' fallback sequence of solution algorithms.

Outputs:
- the SDOF displacement history (relative to the floor) and the spring force;
- the peak |u|;
- the demand at each support: peak displacement Γ·φᵢ·u_peak, and its ratio to the trapeze envelope
  deformations.

The largest support demand, Γ·max φ·u_peak, is the quantity the paper compares with the 3D model.

## 5. 3D verification (paper archetypes)

The paper's full 3D models are stored unchanged in `inputs/models3d/` and include:
- threaded-joint hysteresis;
- Pinching4 trapezes with failure caps;
- both loading directions at once.

They're analysed under a bidirectional floor motion (the two horizontal components of one ground
motion) with the procedure of the `3D_models/*_biron.py` scripts:
- gravity, then mass-proportional damping, 1% at the first (elastic) mode, as in the SDOF;
- Newmark integration in chunks of 100 steps, stopped (collapse) when a braced node exceeds 60 mm;
- the same fallbacks, with each script's own convergence-test types (they differ between models);
- displacements recorded every 0.001 s, as the scripts.

The scripts store a collapsed run as 100 mm and a run that didn't reach the end of the motion as 0.

The reduced-order prediction of each direction (archetype system file, SDOF at the chosen Δc, default
trapezes) is computed under the same component and compared with the paper's measure: **the largest
peak displacement over the braced nodes** of the 3D model.

## Scope and limits

- **Geometry.** One straight main line plus orthogonal branches, one direction per file. No bends,
  elevation changes, loops, or branches off branches. Branches are lumped (no flexibility, no transverse
  braces of their own) and use the main line's pipe section.
- **Trapezes.** All braces of one type are identical.
- **3D verification.** Available only for the 9 archetypes with a model in `inputs/models3d/`; new 3D
  models can't be generated from a system file.
- **Motions.** Only sets described in `motions/motion_sets.yaml`.
- **Run times** (this machine):

  | Analysis | Time |
  |---|---|
  | Pseudo-pushover (50 steps) | 0.1–4 s |
  | One SDOF time history | about 0.3 s |
  | One 3D run, M01–M03 | about 30 s |
  | One 3D run, M29–M31 | about 4 min |

## Validation against the scripts

Full details: [validation.md](validation.md). All comparisons use the engine's default settings.

| Analysis | Compared with | Result | Report |
|---|---|---|---|
| Archetype inputs | the pushover drivers' arguments | identical | [validation/report.md](../validation/report.md) |
| Pseudo-pushover, 18 archetypes | `Pushover2D/pushover_results_*.txt` | all values identical (3 decimals) | [validation/report.md](../validation/report.md) |
| SDOF parameters, 18 archetypes | `Pushover_SDOF/sdof_from_2d.py` | identical (3 decimals) | [validation/report.md](../validation/report.md) |
| SDOF parameters, proposed procedure | `code_proposed_procedure/equivalent_static.py` | identical (rel. < 1e-12) | [validation/report.md](../validation/report.md) |
| SDOF pushover | `Pushover_SDOF/Results*/DispC.out`, `VbaseC.out` | identical (6 digits); M61y's stored files predate its pushover results | [validation.md](validation.md) |
| SDOF time history | `Pushover_SDOF/Results*/*_NLTHA_peak_displacements.csv`, `code_proposed_procedure/Results/NLTHA/` | identical (4 decimals; all 7920 + 150 runs, collapse flags included) | [validation/timehistory_report.md](../validation/timehistory_report.md) |
| 3D verification | `Results/*_biron_DispX/Y.txt` | identical (2 decimals), collapses included | [validation/timehistory_report.md](../validation/timehistory_report.md) |

The SDOF time histories are compared both with the SDOF the scripts analyse and with the engine's own
SDOF (default settings, which read it as the scripts do).
