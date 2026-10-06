"""
Peak displacement and displaced shape of a 3D NLTHA run, stored in the formats of Results/:

  Results/<tag>_DispX.txt, _DispY.txt          peak absolute displacement over the braced nodes (mm)
                                               (rows: record applied in that direction, in Names.txt
                                               order; columns: IM1-IM10; 0 = not available)
  Results/<tag>DispShapeX.npy, DispShapeY.npy  absolute displaced shape at the time step of that peak,
                                               normalized by the reference node (record, node, IM)

Several scripts (or several runs of one script, e.g. different IMs in parallel terminals) can write
to the same files: each update is done under a lock file.
"""
import os
import time
import numpy as np

RESULTS_DIR = os.environ.get("RESULTS3D_DIR",
                             os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Results"))
N_RECORDS = 44
N_IMS     = 10


def peak_and_shape(path, nodes, braced, ref):
    """
    From a Node recorder file ('-time', one displacement column per node of `nodes`): the peak
    absolute displacement over the `braced` nodes, and |u / u_ref| of all nodes at that time step.
    """
    u = np.loadtxt(path, ndmin=2)[:, 1:]
    a = np.abs(u[:, [nodes.index(n) for n in braced]])
    k = int(np.argmax(a.max(axis=1)))
    return float(a[k].max()), np.abs(u[k] / u[k, nodes.index(ref)])


def save_run(tag, direction, row, im, peak, shape, n_nodes):
    """Store one run (row = record index in Names.txt, im = 0-based IM) of one direction."""
    if not 0 <= im < N_IMS:
        raise ValueError(f"IM{im + 1}: the Results files hold IM1-IM{N_IMS}")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    f_peak  = os.path.join(RESULTS_DIR, f"{tag}_Disp{direction}.txt")
    f_shape = os.path.join(RESULTS_DIR, f"{tag}DispShape{direction}.npy")

    with _Lock(os.path.join(RESULTS_DIR, f"{tag}.lock")):
        peaks  = np.loadtxt(f_peak, ndmin=2) if os.path.exists(f_peak) else np.zeros((N_RECORDS, N_IMS))
        shapes = np.load(f_shape) if os.path.exists(f_shape) else np.zeros((N_RECORDS, n_nodes, N_IMS))
        peaks[row, im]     = peak
        shapes[row, :, im] = shape
        _replace(f_peak, lambda f: np.savetxt(f, peaks, fmt="%.2f"))
        _replace(f_shape, lambda f: np.save(f, shapes))
    print(f"Saved {tag} {direction}: record {row + 1}, IM{im + 1}, peak = {peak:.2f} mm")


def _replace(path, write):
    """Write to a temporary file, then swap it in (retrying while another program has the file open)."""
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        write(f)
    for _ in range(100):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            time.sleep(0.1)
    os.replace(tmp, path)


class _Lock:
    """Lock file: created exclusively on entry, removed on exit (stale locks expire after 2 min)."""

    def __init__(self, path):
        self.path = path

    def __enter__(self):
        while True:
            try:
                os.close(os.open(self.path, os.O_CREAT | os.O_EXCL))
                return self
            except FileExistsError:
                try:
                    if time.time() - os.path.getmtime(self.path) > 120:
                        os.remove(self.path)
                        continue
                except FileNotFoundError:
                    continue
                time.sleep(0.05)

    def __exit__(self, *exc):
        os.remove(self.path)
