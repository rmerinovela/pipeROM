import openseespy.opensees as op
import numpy as np
from matplotlib import pyplot as plt
from Functions import*
import os

# Output files (results, layout, figures) next to this script, wherever it is run from
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

#Lumped orthogonal line
# 2D model of the loading along x, consistent with the 3D model (3D_models/M30_biron.py):
#   the analysed line is a branch line through the main, 2D x = 3D y + 5 m (branches from y = -5 to
#   +5 m); hangers at y = -4.5, -1.5, 1.5, 4.5 m; transverse restraints at y = -4.5 (-y branches) and
#   y = 1.5 m (+y branches); 6 branches lumped on the -y side, 9 on the +y side; the main line (32 m) is
#   the orthogonal at the centre, with the longitudinal restraints of the 3D main line.


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
Lpipe  = 10000
npipes = 3
n_mains = 1
n_mains_left  = 6   # mains lumped left of x_center: the 6 branches on the -y side (2D x = 5 m + y)
n_mains_right = 9   # mains lumped right of x_center: the 9 branches on the +y side
x_center      = Lpipe/2   # split location along the main
Dext = 127
Dint = 113
E = 210000
G = 81000.0
rho = 7.85e-9

#Lpipe  = 32000

# Stiff spring coordinates (user input)
#x_stiff_user = np.array([2000.0, 14000.0, 26000.0])

x0_soft = 1000          # not used: explicit hanger positions below
soft_spacing = 3000
x_hangers_user = np.array([500,3500,6500,9500])   # 3D y = -4.5, -1.5, 1.5, 4.5 m

# NEW: stiff hangers defined by mask (0/1)
# Example mask — must be defined by you
stiff_mask = np.array([
    1,0,1,0   # y = -4.5 m (-y branches) and 1.5 m (+y branches)
], dtype=int)


# ------------------------------------------------------------
# Orthogonal branches (generalized)
# ------------------------------------------------------------
x_ortho_user         = np.array([5000])   # the main line (3D y = 0)
L_ortho_user         = np.array([32000])
n_ortho_pipes_user   = np.array([3])
n_ortho_springs_user = np.array([2])   # longitudinal restraints of the 3D main line (8, 26 m)

n_orth = len(x_ortho_user)



# ------------------------------------------------------------
# Initial assumed shape
# ------------------------------------------------------------
# Number of soft hangers
x_soft = x0_soft + soft_spacing * np.arange(len(stiff_mask))
n_hangers = len(x_soft)

# Initial assumed shape (all DOFs = 1)
d0 = np.ones(n_hangers + n_orth)


# ------------------------------------------------------------
# Define pushover displacement steps
# ------------------------------------------------------------
dc_vec = np.linspace(0.1, 50.0, 50)   # example: 20 steps

# Storage list for all results
results = []

# ------------------------------------------------------------
# Initial shape for the FIRST step
# ------------------------------------------------------------
d_init_current = d0.copy()

# ------------------------------------------------------------
# Loop over displacement steps
# ------------------------------------------------------------
for dc in dc_vec:

    print(f"\n--- Running pushover step: Δc = {dc:.3f} ---")
    (
        d_star,
        uy,
        mode1,
        x_nodes_sorted,   # FE x-coordinates sorted
        x_springs,
        x_stiff_out,
        x_d,              # DOF coordinates (unsorted)
        Gamma,
        M_eff,
        Vb,
        mass_ratio,
        order,            # FE sorting index
        f_push,
    ) = iterate_shape_from_static(
        d_init=d0,
        Delta=dc,
        max_iter=50,
        tol=1e-3,
        x0_soft=x0_soft,
        x_hangers_user=x_hangers_user,
        soft_spacing=soft_spacing,
        stiff_mask=stiff_mask,

        # main-line parameters
        Lpipe=Lpipe,
        Dext=Dext,
        Dint=Dint,
        npipes=npipes,
        n_mains=n_mains,
        n_mains_left=n_mains_left,
        n_mains_right=n_mains_right,
        x_center=x_center,
        E=E,
        G=G,
        rho=rho,

        # generalized orthogonal inputs
        x_ortho_user=x_ortho_user,
        L_ortho_user=L_ortho_user,
        n_ortho_pipes_user=n_ortho_pipes_user,
        n_ortho_springs_user=n_ortho_springs_user,
    )
    

    # --------------------------------------------------------
    # Update initial guess for next step
    # --------------------------------------------------------
    d_init_current = d_star.copy()

    # --------------------------------------------------------
    # Scaled displaced shape
    # --------------------------------------------------------
    d_scaled = dc * d_star

    # Value at the maximum-x DOF (last DOF)
    ref = d_scaled[-1]

    # --------------------------------------------------------
    # Normalized displaced shape
    # --------------------------------------------------------
    d_norm = d_scaled / ref

    # --------------------------------------------------------
    # Equivalent SDOF displacement
    # --------------------------------------------------------
    u_sdof = dc / (Gamma*np.max(d_norm))

    # --------------------------------------------------------
    # Build one row of results:
    # [dc, Gamma, M_eff, Vb, mass_ratio, u_sdof, d_norm..., d_scaled...]
    # --------------------------------------------------------
    row = (
        [dc, Gamma, M_eff, Vb, mass_ratio, u_sdof]
        + list(d_norm)
        + list(d_scaled)
        + list(f_push)
    )

    results.append(row)


# ------------------------------------------------------------
# Convert to array and write to txt file (no headers)
# ------------------------------------------------------------
results_arr = np.array(results, float)

np.savetxt(
    os.path.join(SCRIPT_DIR, "pushover_results_M30x.txt"),
    results_arr,
    fmt="%.3f",
)

# Layout of the 2D model, used to build the SDOF (Pushover_SDOF)
save_pushover_layout(
    os.path.join(SCRIPT_DIR, "pushover_layout_M30x.json"),
    x_d,
    x_stiff_out,
    n_ortho_springs_user,
    n_mains_left,
    n_mains_right,
    x_center,
)

# ------------------------------------------------------------
# Displaced shape at critical displacement vs NLTHA mean shape
# ------------------------------------------------------------
# NLTHA points on the analysed main: (label, point indices, 2D x [mm])
nltha_lines = [
    ('(line at x = 7.05 m)',
     [55, 56, 57, 7, 31, 32, 33],
     [0, 500, 3500, 5000, 6500, 9500, 10000]),
    ('(line at x = 11.55 m)',
     [58, 59, 60, 10, 34, 35, 36],
     [0, 500, 3500, 5000, 6500, 9500, 10000]),
]

plot_shape_vs_nltha(
    results_arr,
    x_d,
    model="M30",
    direction="X",
    nltha_lines=nltha_lines,
    dc_target=12.0,
    level=-1,                     # strongest NLTHA intensity
    x_supports=x_stiff_out,
    out_png=os.path.join(SCRIPT_DIR, "shape_comparison_M30x.png"),
)
