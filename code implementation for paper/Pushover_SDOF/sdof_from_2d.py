import os
import json
import numpy as np

PUSHOVER2D_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Pushover2D")


def load_sdof_params(tag, dc_target=12.0, pushover_dir=PUSHOVER2D_DIR):
    """
    Equivalent SDOF parameters from the 2D pushover of archetype `tag` (e.g. "M01x").

    Reads Pushover2D/pushover_results_<tag>.txt (rows: [dc, Gamma, M_eff, Vb,
    mass_ratio, u_sdof, d_norm..., d_scaled..., f_push...]) at the step closest
    to dc_target, and Pushover2D/pushover_layout_<tag>.json.

    Returns a dict with
      mass   : M_eff (t)
      Gamma  : participation factor
      dc     : pushover step actually used (mm)
      trans  : list of (phi, scale) for the transverse springs (stiff hangers)
      long   : list of (phi, scale) for the longitudinal springs (orthogonals)
    phi is the normalized displacement (d_norm) at the spring's DOF; scale is
    the number of identical springs it represents (mains lumped on that side
    of x_center for the hangers, alpha * n springs for the orthogonals), so
    that the SDOF base shear matches the 2D Vb.
    """
    results = np.loadtxt(os.path.join(pushover_dir, f"pushover_results_{tag}.txt"), ndmin=2)
    with open(os.path.join(pushover_dir, f"pushover_layout_{tag}.json")) as f:
        layout = json.load(f)

    x_d   = np.asarray(layout["x_d"], float)
    n_dof = len(x_d)
    if results.shape[1] != 6 + 3 * n_dof:
        raise ValueError(f"{tag}: results have {results.shape[1]} columns, "
                         f"layout has {n_dof} DOFs. Rerun the 2D pushover.")

    i_step = np.argmin(np.abs(results[:, 0] - dc_target))
    row    = results[i_step]

    # Back to DOF order (some drivers save the shape sorted by x)
    d_norm = np.empty(n_dof)
    d_norm[layout["column_order"]] = row[6:6 + n_dof]

    # Mains lumped at x (same rule as build_model's n_mains_at)
    nL, nR, xc = layout["n_mains_left"], layout["n_mains_right"], layout["x_center"]

    def n_mains_at(x, tol=1e-3):
        if x < xc - tol:
            return float(nL)
        if x > xc + tol:
            return float(nR)
        return 0.5 * (nL + nR)

    trans = [(d_norm[i], n_mains_at(x_d[i])) for i in layout["i_stiff"]]
    long_ = [(d_norm[i], layout["alpha"] * ns)
             for i, ns in zip(layout["i_ortho"], layout["n_ortho_springs"])]

    for phi, _ in trans + long_:
        if phi <= 1e-3:
            raise ValueError(f"{tag}: non-positive shape value {phi:.3g} at a spring DOF "
                             f"(dc = {row[0]:.2f} mm); choose another dc_target.")

    params = {
        "mass":  row[2],
        "Gamma": row[1],
        "dc":    row[0],
        "trans": trans,
        "long":  long_,
    }

    print(f"{tag} SDOF from 2D pushover at dc = {row[0]:.2f} mm: "
          f"M_eff = {params['mass']:.3f}, Gamma = {params['Gamma']:.3f}")
    print("  transverse (phi, n):  ", [(round(p, 3), s) for p, s in trans])
    print("  longitudinal (phi, n):", [(round(p, 3), s) for p, s in long_])
    return params


def plot_pushover_comparison(tag, disp_file, force_file, dc_target=12.0, cyclic=False,
                             out_png=None, pushover_dir=PUSHOVER2D_DIR):
    """
    SDOF pushover (OpenSees recorders) vs the 2D adaptive pushover of Pushover2D.

    disp_file  : SDOF displacement recorder (DispC*.out)
    force_file : base reaction recorder (VbaseC*.out); force = -reaction
    The 2D curve is Vb vs u_sdof = u_ref / Gamma, the displacement of the SDOF,
    mirrored for negative displacements when cyclic=True. The marker shows the
    2D step the SDOF was built from (dc_target).
    """
    import matplotlib.pyplot as plt

    D = np.loadtxt(disp_file).ravel()
    F = -np.loadtxt(force_file).ravel()
    n = min(len(D), len(F))

    results = np.loadtxt(os.path.join(pushover_dir, f"pushover_results_{tag}.txt"), ndmin=2)
    u2d = np.concatenate([[0.0], results[:, 5]])
    V2d = np.concatenate([[0.0], results[:, 3]])
    i_ref = np.argmin(np.abs(results[:, 0] - dc_target))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(D[:n], F[:n] / 1000, color='#365c8d', lw=1.2,
            label='SDOF ' + ('cyclic pushover' if cyclic else 'pushover'))
    ax.plot(u2d, V2d / 1000, 'k--', lw=2, label='2D adaptive pushover')
    if cyclic:
        ax.plot(-u2d, -V2d / 1000, 'k--', lw=2)
    ax.plot(results[i_ref, 5], results[i_ref, 3] / 1000, 'o', color='#d1495b', ms=8,
            label=rf'SDOF calibration ($\Delta_c$ = {results[i_ref, 0]:.1f} mm)')

    ax.axhline(0, color='0.7', lw=0.8)
    ax.axvline(0, color='0.7', lw=0.8)
    ax.set_xlabel('SDOF displacement (mm)', fontsize=13)
    ax.set_ylabel('Base shear (kN)', fontsize=13)
    ax.set_title(f'{tag}: SDOF vs 2D pushover')
    ax.legend(fontsize=10)
    fig.tight_layout()

    if out_png is not None:
        fig.savefig(out_png, dpi=200)
        print(f"Saved pushover comparison: {out_png}")
    return fig, ax
