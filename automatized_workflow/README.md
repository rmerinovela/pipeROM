# Automatized workflow

The workflow the app follows when the user describes a piping system that is not one of the paper's
archetypes: all results are created from scratch, with the same steps used for the paper.

Status: **in progress**. Step 1 (3D model) is done.

## Why one description of the whole system

For the paper, the 3D model of each archetype and its two directional reduced-order models were written
by hand. They disagree in several archetypes, e.g. the transverse restraints of M30y and of the old M62y
are one hanger apart in the 3D model and in the 2D stiff mask. Here everything is derived from a single
description of the system in plan, so the 3D model and the reduced-order models agree by construction.

## Steps

| # | Step | Produces | Built on | Status |
|---|---|---|---|---|
| 0 | System description | one YAML file: main line, hangers, branches (position, side, length), pipe, restraints per line and direction, trapezes | format below | done |
| 1 | 3D model | OpenSees model in the format of `inputs/models3d/*.json` (commands + `nodesX/Y`, `nodesXt/Yt`, `nodesXp/Yp`) | `model3d.py` | done |
| 2 | Reduced-order system files | `<name>x.yaml`, `<name>y.yaml` (format of `inputs/archetypes/`), plus the map from the 2D line to the 3D shape nodes (`nltha_lines`) | new | to do |
| 3 | Pseudo-pushover and SDOF | pushover and SDOF of both directions | `piperom.pushover`, `piperom.sdof` | to do |
| 4 | SDOF time histories | peak SDOF displacements, IM1–IM10 | `piperom.timehistory` | to do |
| 5 | 3D time histories | peak displacement at the restrained nodes and displaced shape at the peak, per record and IM, in parallel | `piperom.verification3d`; conventions of `code implementation for paper/3D_models/results3d.py`, `run_3d.py` | to do |
| 6 | Comparison | RE vs intensity, record-by-record peaks, displaced shapes (static vs NLTHA) | `code implementation for paper/Results/nltha_comparison.py` | to do |
| 7 | Driver | runs steps 1–6 from the description; called by the app | new | to do |

## System description

```yaml
name: MySystem                               # units N, mm, s
description: free text                       # optional
pipe: {...}                                  # optional, as in the system files (defaults: inputs/defaults/system.yaml)
main_line:                                   # along x, from x = 0
  length: 36000
  n_pipes: 3
  hangers: {first: 1000, spacing: 3000}      # or {positions: [...]}; optional end_clearance
  braces:                                    # restraints, as in the app: {mask: [...]}, {positions: [...]} or {count: n}
    transverse: {count: 4}                   #   along y
    longitudinal: {count: 3}                 #   along x
branches:                                    # along y, framing into the main line
  - x: 36000                                 # tee position (0 < x <= length)
    side: +y                                 # +y or -y
    length: 18000
    n_pipes: 3                               # optional, default: main line's
    hangers: {first: 1000, spacing: 3000}    # optional, default: main line's; distances from the tee
    braces:
      transverse: {count: 2}                 #   along x
      longitudinal: {count: 2}               #   along y
trapezes: {transverse: default, longitudinal: default}   # optional, as in the system files
```

Examples: `examples/M02.yaml` (counts), `examples/M62.yaml` (positions of the repo's M62 3D model).

## Step 1: 3D model (`model3d.py`)

```bash
python automatized_workflow/model3d.py automatized_workflow/examples/M02.yaml -o M02_auto.json
```

or `build_model3d(load_description(path))` from Python. The model is the one of the generator
`SuspendedNSEs/GeneralModel/M02.py` (the paper's 3D models, generalized), with the fixes of 2026-10-05:
branch joints measured from the tee and none at a free end; one tee node per branch position (an
existing main-line node at that position is the tee); gravity from the node masses; node lists for the
NLTHA output. Given the same restraint positions, the two give identical OpenSees commands and node lists
(checked on M02, M62 and a three-branch layout with a tee at a hanger and branches on both sides); counts
are placed with a different rule (below).

Choices (2026-10-05):
- **Restraints** as in the app: mask, positions or count.
- **Placement of a count**: the app's rule (`piperom.inputs.evenly_spaced_brace_mask`, the
  `compute_stiff_mask` of the paper's `equivalent_static.py`): spacing L/n, starting half a spacing from
  the start of the line, each snapped to the nearest free hanger. The same rule for transverse and
  longitudinal restraints and for every line (branches measured from the tee). Counts are resolved here
  once and step 2 will pass the resulting positions to the reduced-order files. (The generator used
  spacing (L − x0)/n from the first hanger x0, which picks different hangers in about half the cases;
  replaced 2026-10-05.)
- **Hangers**: the piperom grid (first, spacing, last at least `end_clearance` from the end), so the 3D
  and reduced-order models have the same hangers.
- **Pipe section** from the diameters (`piperom.static_model._section`). The paper's 3D models use typed
  values (A = 7761.3 mm², I = 17.86·10⁶ mm⁴), 25% stiffer in bending.
- **Mass**: (pipe + fluid) per unit length × mass factor, as in piperom (default mass factor 1.35, fluid
  density = density / 7.8), shared equally by the hangers of each line.
- **Trapezes**: the piperom trapeze files (the same Pinching4 parameters as the paper's 3D models).
- Branch elements use the branch's number of pipes (the generator used the main line's for all lines).

With the paper's restraint positions and section properties the model reproduces the repo's 3D models of
M02 and M62 (same nodes, elements, joints, hangers and node lists; periods within 0.4%).

## Open questions

- **Reduced-order model along the main line (step 2).** For loading along x, the paper's 2D model
  analyses the branches lumped into one line (M62x: the two 52 m branches as one line of 6 pipes, the main
  line as the orthogonal branch; M29x: one 12 m line of 18 pipes for its 15 short branches). The rule for
  general layouts (branches of different lengths or numbers of pipes, both sides of the main line) is
  needed.
- **Layouts**: one straight main line with branches perpendicular to it; uniform hanger spacing or
  explicit positions on each line. Irregular grids such as M29–M31 should be expressible as hanger
  positions (main line shifted to start at x = 0); not tested yet.
