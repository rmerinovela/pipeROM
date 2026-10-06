import os
import csv

COLUMNS = ["IM", "record", "peak_disp_mm", "completed", "collapsed"]


def save_peak_displacements(path, im, records, peak_disp, completed, collapsed=None):
    """
    Write the peak absolute SDOF displacements of one IM level to the
    archetype's results file (one row per record), replacing any previous
    rows of that IM. Rows of the other IMs are kept. collapsed: True where a
    braced node exceeded the collapse displacement and the analysis was stopped
    (files written before this column read as 0).
    """
    if collapsed is None:
        collapsed = [False] * len(records)
    rows = []
    if os.path.exists(path):
        with open(path, newline="") as f:
            rows = [{**r, "collapsed": r.get("collapsed") or 0}
                    for r in csv.DictReader(f) if int(r["IM"]) != int(im)]

    rows += [{"IM": int(im), "record": str(rec), "peak_disp_mm": f"{d:.4f}",
              "completed": int(bool(c)), "collapsed": int(bool(k))}
             for rec, d, c, k in zip(records, peak_disp, completed, collapsed)]
    rows.sort(key=lambda r: int(r["IM"]))     # stable: record order kept within an IM

    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)

    n_coll = sum(bool(k) for k in collapsed)
    n_bad  = sum(not c and not k for c, k in zip(completed, collapsed))
    print(f"Saved {len(records)} peak displacements of IM{im} to {os.path.basename(path)}"
          + (f" ({n_coll} collapsed)" if n_coll else "")
          + (f" ({n_bad} records did not reach the end)" if n_bad else ""))
