"""
Static pushover of the 3D models, to compare with the 2D adaptive pushover and the equivalent SDOF.

The model is the one of the NLTHA (<M>_biron.py), built from inputs/models3d/<M>.json (the same OpenSees
commands, recorded from the script). After gravity, a lateral load proportional to the nodal masses
(uniform acceleration) is applied in X or Y and increased under displacement control at the control
node: the 3D node where the 2D model of that direction is normalized (its reference DOF, as in the
displaced-shape plots of Results/nltha_comparison.py). Base shear = applied lateral load (load factor x
total mass), equal to the support reactions in a static analysis (the reactions of the fixed nodes aren't
available directly: the trapezes reach them through rigid links).

    python pushover_3d.py M01                 both directions, M01
    python pushover_3d.py M01 M29 --dir X     direction X of M01 and M29
    python pushover_3d.py all                 every model, in parallel

Output: Results/<M>_biron_pushover<X|Y>.txt, columns: control displacement (mm), base shear (N).
"""
import os
import sys
import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import openseespy.opensees as op

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
REPO_DIR    = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
RESULTS_DIR = os.path.join(SCRIPT_DIR, '..', 'Results')
sys.path.insert(0, REPO_DIR)
sys.path.insert(0, RESULTS_DIR)

MODELS = ['M01', 'M02', 'M03', 'M29', 'M30', 'M31', 'M61', 'M62', 'M63']
TARGET = 60.0     # control displacement (mm)
DSTEP  = 0.05     # displacement increment (mm)


def control_node(model, direction):
    """3D node at the reference DOF of the 2D model of this direction (nltha_comparison.static_shape)."""
    from nltha_comparison import load_geometry, static_shape
    geom = load_geometry(model)
    return geom['nodes'][static_shape(model, direction, geom)['i_ref']]


def pushover(model, direction, target=TARGET, dstep=DSTEP):
    from piperom.verification3d import build_model3d, load_model3d
    m = load_model3d(model)
    ctrl, dof = control_node(model, direction), 1 if direction == 'X' else 2

    build_model3d(m)
    # Gravity, as in the NLTHA scripts
    op.constraints('Transformation')
    op.numberer('RCM')
    op.system('BandGeneral')
    op.test('NormDispIncr', 1.0e-8, 200)
    op.algorithm('Newton')
    op.integrator('LoadControl', 0.1)
    op.analysis('Static')
    if op.analyze(10) != 0:
        raise RuntimeError(f'{model}: gravity analysis failed')
    op.loadConst('-time', 0.0)

    # Lateral load proportional to the nodal masses in the loading direction
    op.timeSeries('Linear', 2)
    op.pattern('Plain', 200, 2)
    total_mass = 0.0
    for nd in op.getNodeTags():
        mass = op.nodeMass(nd, dof)
        if mass > 0:
            load = [0.0] * 6
            load[dof - 1] = mass
            op.load(nd, *load)
            total_mass += mass

    op.wipeAnalysis()
    op.constraints('Transformation')
    op.numberer('RCM')
    op.system('UmfPack')
    op.test('EnergyIncr', 1.e-4, 1000)
    op.algorithm('Newton')
    op.integrator('DisplacementControl', ctrl, dof, dstep)
    op.analysis('Static')

    u, vb = [0.0], [0.0]
    for _ in range(int(round(target / dstep))):
        ok = op.analyze(1)
        for algorithm, n_sub in ((('Newton', '-initial'), 2), (('Broyden', 50), 2), (('NewtonLineSearch',), 10)):
            if ok == 0:
                break
            op.test('EnergyIncr', 1.e-3, 5000)
            op.algorithm(*algorithm)
            op.integrator('DisplacementControl', ctrl, dof, dstep / n_sub)
            ok = op.analyze(n_sub)
            op.integrator('DisplacementControl', ctrl, dof, dstep)
            op.test('EnergyIncr', 1.e-4, 1000)
            op.algorithm('Newton')
        if ok != 0:
            print(f'{model} {direction}: no convergence at {op.nodeDisp(ctrl, dof):.2f} mm, stopped')
            break
        u.append(op.nodeDisp(ctrl, dof))
        vb.append(op.getLoadFactor(200) * total_mass)
    op.wipe()

    path = os.path.join(RESULTS_DIR, f'{model}_biron_pushover{direction}.txt')
    np.savetxt(path, np.column_stack([u, vb]), fmt='%.6g',
               header=f'{model} 3D pushover in {direction}: mass-proportional load, control node {ctrl}\n'
                      f'control displacement (mm), base shear (N)')
    print(f'Saved {os.path.basename(path)}: {len(u) - 1} steps, Vb at {u[-1]:.1f} mm = {vb[-1] / 1e3:.1f} kN')
    return path


def _run(args):
    return pushover(*args)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('models', nargs='+', help="model names (e.g. M01) or 'all'")
    ap.add_argument('--dir', choices=['X', 'Y'], help='one direction only (default: both)')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = ap.parse_args()
    models = MODELS if args.models == ['all'] else args.models
    tasks = [(m, d) for m in models for d in ([args.dir] if args.dir else ['X', 'Y'])]
    with ProcessPoolExecutor(min(args.workers, len(tasks))) as pool:
        list(pool.map(_run, tasks))


if __name__ == '__main__':
    main()
