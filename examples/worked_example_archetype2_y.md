# Worked example: archetype 2, y direction

This example takes archetype 2 of the paper, loaded in the y direction, through Phases 1 and 2 of the
reduced-order procedure (Sections 3.1 and 3.2 of the paper), with the actual numbers of each step. Step and
equation numbers refer to the paper. Every number can be reproduced with the notebook
[`worked_example_archetype2_y.ipynb`](worked_example_archetype2_y.ipynb), which computes each step explicitly
and checks it against `piperom`. Units are N, mm and t.

## 1. Input data

The main line is 36 m long and consists of three parallel 127 mm water-filled pipes. Gravity hangers are spaced
at 3 m, starting 1 m from its left end. In the y direction the main line is the transverse pipeline, with
transverse restraints at the hangers at x = 1, 13, 22 and 34 m. The 18 m branch that frames into its right end
(x = 36 m) is the longitudinal pipeline, with two longitudinal restraints (system file:
[`inputs/archetypes/M02y.yaml`](../inputs/archetypes/M02y.yaml)).

The mass of pipes and water (62.1 + 30.3 = 92.4 kg/m) is increased by 35% to account for fittings and the
weight of the trapezes, giving 124.8 kg/m. Table 1 gives the restraint backbones used by the equivalent static
analysis: the first-yield and second points of the restraints' envelopes
([`inputs/trapezes/`](../inputs/trapezes/)), with a post-yield slope of 1% of the elastic one.

**Table 1.** Restraint backbones.

| Restraint | $d_1$ (mm) | $F_1$ (N) | $d_2$ (mm) | $F_2$ (N) | $k_1$ (N/mm) |
|---|---|---|---|---|---|
| Transverse | 10 | 6000 | 17 | 9000 | 600 |
| Longitudinal | 12 | 7500 | 24 | 10000 | 625 |

## 2. Phase 1: equivalent single-pipeline model

**Step 1: lump the longitudinal line.** The branch becomes a mass $m_{L,1}$ = 124.8 kg/m × 18 m = 2.246 t on a
spring $k_{L,1}$, its two longitudinal restraints in parallel: (12 mm, 15 000 N), (24 mm, 20 000 N).

**Step 2.1: segments.** The branch frames into the end of the transverse line, so the line has one segment
($t$ = 1, from 0 to 36 m).

**Step 2.2: masses.** The masses are lumped at the hangers and at the junction using tributary lengths
(Table 2). The tributary mass carried by the longitudinal line is $m_{T|L,1}$ = 0.125 t, from 35 m (half-way to
the restraint at 34 m) to 36 m. The first 0.5 m at the free left end (0.06 t, under 1% of the mass) lumps onto the
end node, which is not a point of support, and is not part of the procedure.

**Table 2.** Lumped masses of the transverse line.

| x (m) | 1 | 4 to 31 (10 hangers) | 34 | 36 (junction) |
|---|---|---|---|---|
| Tributary length (m) | 2.0 | 3.0 each | 2.5 | 1.0 |
| Mass (t) | 0.250 ($m_{T,1,1}$) | 0.374 each | 0.312 | 0.125 ($m_{T\|L,1}$) |

**Step 2.3: parallel transverse lines.** Not needed: there is only one transverse line in this direction.

**Step 3: equivalent linear static analysis.** The target displacement $\Delta_c$ = 12 mm is assigned to the
largest displacement of the transverse line; with it, the most-displaced transverse restraint reaches
approximately its first-yield displacement (10 mm). The initial guess is a uniform shape (step a: $\Delta$ = 12 mm
at all points, with $\Delta_{L,1} = \Delta_{T|L,1}$). In the first iteration:

- **b)** Each transverse restraint at 12 mm carries 6857 N (secant stiffness 571 N/mm), so $V_{T,tot,1}$ = 27 429 N.
  The longitudinal spring carries $V_{L,tot,1}$ = 15 000 N.
- **c)** Eq. 3: $V_{L,1}$ = 2.246 / (2.246 + 0.125) × 15 000 = 14 210 N stays on the longitudinal line.
  Eq. 4 ($\gamma$ = 1, since the junction is at the end of segment 1): $V_{T|L,1,1}$ = 790 N is passed to the
  transverse line.
- **d)** Eq. 6: $V_{tot,1}$ = 27 429 + 790 = 28 218 N, distributed in proportion to $m\Delta$ (Eqs. 7 and 8). With a
  uniform $\Delta$ this is proportional to the masses: 2385 N at each 3 m hanger.
- **e)** The linear static analysis (elastic beam on the secant restraint springs, plus the lumped longitudinal
  line loaded by $V_{L,1}$) gives $\Delta'$, with an RMS change of 0.167 (normalised shape) from the initial guess.

The procedure converges after 9 iterations (RMS change of the normalised shape below 0.001, i.e. 0.012 mm).
Table 3 lists the converged shape. The largest displacement, 12 mm, is at the unbraced hanger at 28 m. The
most-displaced transverse restraint (22 m) reaches 10.0 mm, its first yield; all other restraints remain elastic.

**Table 3.** Converged displaced shape and restraint forces.

| x (m) | Point | $\Delta$ (mm) | Restraint force (N) |
|---|---|---|---|
| 1 | transverse restraint | 7.48 | 4489 |
| 4 | hanger | 10.22 | |
| 7 | hanger | 11.62 | |
| 10 | hanger | 11.49 | |
| 13 | transverse restraint | 9.87 | 5922 |
| 16 | hanger | 11.41 | |
| 19 | hanger | 11.46 | |
| 22 | transverse restraint | 10.01 | 6005 |
| 25 | hanger | 11.78 | |
| 28 | hanger | 12.00 | |
| 31 | hanger | 10.64 | |
| 34 | transverse restraint | 7.90 | 4737 |
| 36 | junction (longitudinal line) | 7.84 | 9802 |

## 3. Phase 2: equivalent SDOF system

**Step 1: effective mass and equivalent displacement.** With $\mathbf{m}_{full}$ and $\boldsymbol{\Delta}_{full}$,
which include $m_{L,1}$ = 2.246 t at $\Delta_{L,1} = \Delta_{T|L,1}$ = 7.84 mm,
$\sum m_l \Delta_l$ = 64.28 t·mm and $\sum m_l \Delta_l^2$ = 638.3 t·mm². Eqs. 9 and 10 give

$$m_{eff} = \frac{(\sum m_l \Delta_l)^2}{\sum m_l \Delta_l^2} = \frac{64.28^2}{638.3} = 6.474\ \text{t}\quad(97\%\ \text{of the total mass of 6.676 t}),
\qquad \Delta_D = \frac{\sum m_l \Delta_l^2}{\sum m_l \Delta_l} = \frac{638.3}{64.28} = 9.93\ \text{mm}.$$

**Step 2: non-linear model of the SDOF.** Each restraint's hysteretic model is moved to SDOF coordinates with
Eq. 11, $\tilde V_l(\Delta) = V_l(\Delta\,\Delta_l/\Delta_D)$: its envelope deformations are divided by
$\Delta_l/\Delta_D$ and its forces are unchanged (Table 4); the cyclic parameters are ratios and stay as they are.
The SDOF spring is the sum of the five transformed springs in parallel.

**Table 4.** Restraints in SDOF coordinates (first-yield and second envelope points).

| Restraint | $\Delta_l$ (mm) | $\Delta_l/\Delta_D$ | $d$ (mm) at 6000 N / 15 000 N | $d$ (mm) at 9000 N / 20 000 N |
|---|---|---|---|---|
| Transverse, 1 m | 7.48 | 0.754 | 13.27 | 22.56 |
| Transverse, 13 m | 9.87 | 0.994 | 10.06 | 17.10 |
| Transverse, 22 m | 10.01 | 1.008 | 9.92 | 16.86 |
| Transverse, 34 m | 7.90 | 0.795 | 12.58 | 21.38 |
| Longitudinal line (2 restraints) | 7.84 | 0.790 | 15.20 | 30.39 |

The resulting equivalent SDOF system ($m_{eff}$ = 6.474 t, $\Delta_D$ = 9.93 mm and the spring of Table 4) is the
one compared with the MDOF model of archetype 2 in the paper (force-displacement curves, displaced shape in the y
direction and fragility functions). The SDOF used there was derived at the pushover step closest to 12 mm
($\Delta_c$ = 12.32 mm); its $m_{eff}$ and $\Delta_l/\Delta_D$ agree with the values above to within 0.3%.

In the SDOF files written by `piperom sdof`, the same transformation is expressed with the participation factor
$\Gamma$ = 0.790 and the shape $\phi$ normalised at the longitudinal line: $\Delta_l/\Delta_D = \Gamma\,\phi_l$.
