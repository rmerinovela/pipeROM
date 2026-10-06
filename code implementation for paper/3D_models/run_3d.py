"""
Run the NLTHA of a 3D model in parallel: one process per (intensity level, ground motion), each
running both orientations of that ground motion. Results go to Results/ as in the single scripts.

    python run_3d.py M62_biron_mask2D                       all IMs (1-10) and ground motions (1-22)
    python run_3d.py M62_biron --ims 1-5 --workers 16       IM1-IM5 on 16 cores
    python run_3d.py M29_biron --ims 6 --gms 1-4 8
    python run_3d.py M62_biron_mask2D --skip-done           only runs missing from Results/ (resume)

The output of a task is kept in logs/ only if it fails.
"""
import os
import re
import sys
import time
import argparse
import subprocess
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from results3d import RESULTS_DIR  # noqa: E402


def numbers(tokens):
    """['1-3', '7'] -> [1, 2, 3, 7]"""
    out = []
    for t in tokens:
        a, _, b = t.partition('-')
        out += list(range(int(a), int(b or a) + 1))
    return sorted(set(out))


def results_tag(script):
    with open(script, encoding='utf-8') as f:
        return re.search(r"^RESULTS_TAG = '(\w+)'", f.read(), re.M).group(1)


def done(tag, im, gm):
    """Both records of the ground motion stored (non-zero) in both directions at this IM."""
    for d in 'XY':
        path = os.path.join(RESULTS_DIR, f'{tag}_Disp{d}.txt')
        if not os.path.exists(path):
            return False
        peaks = np.loadtxt(path, ndmin=2)
        if not (peaks[2*gm - 2, im - 1] > 0 and peaks[2*gm - 1, im - 1] > 0):
            return False
    return True


def run(script, im, gm):
    t0 = time.time()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    p = subprocess.run([sys.executable, script, str(im), '--gm', str(gm)], cwd=SCRIPT_DIR, env=env,
                       capture_output=True, text=True)
    saved = re.findall(r'^Saved .*$', p.stdout, re.M)
    ok = p.returncode == 0 and len(saved) == 4
    if not ok:
        os.makedirs(os.path.join(SCRIPT_DIR, 'logs'), exist_ok=True)
        log = os.path.join(SCRIPT_DIR, 'logs', f'{os.path.basename(script)[:-3]}_IM{im}_GM{gm}.log')
        with open(log, 'w') as f:
            f.write(p.stdout + '\n' + p.stderr)
    return ok, time.time() - t0, saved


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('model', help='script name in 3D_models, e.g. M62_biron or M62_biron_mask2D')
    ap.add_argument('--ims', nargs='+', default=['1-10'], help='intensity levels, e.g. 1-5 8 (default 1-10)')
    ap.add_argument('--gms', nargs='+', default=['1-22'], help='ground motions, e.g. 1-4 (default 1-22)')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument('--skip-done', action='store_true', help='skip runs already stored in Results/')
    args = ap.parse_args()

    script = os.path.join(SCRIPT_DIR, args.model if args.model.endswith('.py') else args.model + '.py')
    tag = results_tag(script)
    # Strongest IMs first: they take longest, so the pool finishes evenly
    tasks = [(im, gm) for im in sorted(numbers(args.ims), reverse=True) for gm in numbers(args.gms)]
    if args.skip_done:
        tasks = [t for t in tasks if not done(tag, *t)]
    print(f'{os.path.basename(script)}: {len(tasks)} tasks (IM x ground motion) on {args.workers} workers '
          f'-> Results/{tag}_*', flush=True)

    t0, failed = time.time(), []
    with ThreadPoolExecutor(args.workers) as pool:
        futures = {pool.submit(run, script, im, gm): (im, gm) for im, gm in tasks}
        for k, fut in enumerate(as_completed(futures), 1):
            im, gm = futures[fut]
            ok, dt, saved = fut.result()
            peaks = ', '.join(re.sub(r'^Saved \S+ (\w): record (\d+), IM\d+, peak = ', r'\1 r\2: ', s) for s in saved)
            status = 'ok' if ok else 'FAILED (see logs/)'
            note = ' [a record did not reach the end of the motion: stored as 0]' if any(s.endswith('peak = 0.00 mm') for s in saved) else ''
            if any(s.endswith('peak = 100.00 mm') for s in saved):
                note += ' [collapse: restrained node above 60 mm, stored as 100 mm]'
            print(f'[{k}/{len(tasks)}] IM{im} GM{gm} {status} in {dt:.0f} s  {peaks}{note}', flush=True)
            if not ok:
                failed.append((im, gm))

    print(f'Done in {(time.time() - t0) / 60:.1f} min' +
          (f'; failed: {failed} (rerun with --skip-done)' if failed else ''))


if __name__ == '__main__':
    main()
