import openseespy.opensees as op
import numpy as np
from matplotlib import pyplot as plt
import sys
import os

# The 2D model (build_model, iterate_shape_from_static, save_pushover_layout) is the one of the paper's
# pushover analyses, Pushover2D/Functions.py, so both use the same procedure
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, '..', 'Pushover2D'))
from Functions import iterate_shape_from_static, save_pushover_layout

def compute_stiff_mask(Lpipe, x0, soft_spacing, nT):
    """
    Place nT supports with spacing ~ Lpipe/nT, starting near Lpipe/(2*nT),
    then snap each to nearest unused hanger.
    """

    # 1. Generate hanger coordinates, as build_model (Pushover2D/Functions.py): from x0 every
    #    soft_spacing, the last one at least 1 m from the right end
    Nmax = int((Lpipe - 1000.0 - x0) // soft_spacing)
    hangers = x0 + soft_spacing * np.arange(Nmax + 1, dtype=float)
    hangers = hangers[hangers <= Lpipe - 1000.0 + 1e-6]
    nhangers = len(hangers)

    # 2. Target support positions based on Lpipe/nT
    spacing = Lpipe / nT
    first = spacing / 2.0
    ideal_supports = first + spacing * np.arange(nT)

    # 3. Snap each ideal support to nearest unused hanger
    snapped_indices = []
    used = set()

    for x in ideal_supports:
        d = np.abs(hangers - x)
        order = np.argsort(d)
        for idx in order:
            if idx not in used:
                snapped_indices.append(idx)
                used.add(idx)
                break

    snapped_indices = np.array(snapped_indices, dtype=int)

    # 4. Build stiff_mask
    stiff_mask = np.zeros(nhangers, dtype=int)
    stiff_mask[snapped_indices] = 1

    return stiff_mask, hangers, snapped_indices



# Driver
# ------------------------------------------------------------
# Input Parameters
# ------------------------------------------------------------
Lpipe  = 36000 # length of transverse pipeline
npipes = 3     # number of pipes in transverse pipeline 
n_mains = 1    # number of parallel mains lumped into the analysed one
n_mains_left  = n_mains   # mains lumped left of x_center
n_mains_right = n_mains   # mains lumped right of x_center
x_center      = Lpipe/2   # split location along the main
Dext = 127     # external diameter of pipe
Dint = 113     # internal diameter of pipe
E = 210000     # modulus of elasticity of pipe material
G = 81000.0    # shear modulus of pipe material
rho = 7.85e-9  # density of pipe material
nT = 4         # number of transverse supports on transverse pipeline
nL = 3         # number of longitudinal supports on longitudinal pipelines

x0_soft =1000       # first distance from left end of transverse line to first gravity hanger
soft_spacing = 3000 # spacing of gravity hangers (seismic supports are always at locations of gravity hangers)

# stiff hangers defined by mask (0/1)
stiff_mask, hangers, snapped_indices = compute_stiff_mask(
    Lpipe, x0_soft, soft_spacing, nT
)


# ------------------------------------------------------------
# Orthogonal branches (generalized)
# ------------------------------------------------------------
# the number of entries all arrays corresponds to the total number of longitudinal lines framing into the transverse one -- the length of all the following arrays must be equal
x_ortho_user         = np.array([36000])   # distance from left end of transverse pipeline where longitudindal lines are located
L_ortho_user         = np.array([18000])   # length of longitudinal pipelines
n_ortho_pipes_user   = np.array([3])       # numbers of pipes in each longitudinal pipeline
n_ortho_springs_user = np.array([nL])      # number of longitudinal supports on each longitudinal pipeline

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
# Define target displacement of equivalent static procedure
# ------------------------------------------------------------
dc = 12   # target displacement

# Storage list for all results
results = []

# ------------------------------------------------------------
# Initial shape for the FIRST step
# ------------------------------------------------------------
d_init_current = d0.copy()

# ------------------------------------------------------------
# Loop over displacement steps
# ------------------------------------------------------------

print(f"\n--- Running Equivalent Static Procedure: Δc = {dc:.1f} ---")
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
# Save to file
# --------------------------------------------------------


# Next to this script, wherever it is run from; read by NLTHA_SDOF.py
outdir = os.path.join(SCRIPT_DIR, "Results", "EquivStatic")
os.makedirs(outdir, exist_ok=True)
tag = f"equiv_static_nT{nT}_nL{nL}"

# Same files as the paper's pushover analyses (Pushover2D/pushover_results_<tag>.txt and
# pushover_layout_<tag>.json), here with the single step dc, so that the SDOF is built with
# Pushover_SDOF/sdof_from_2d.py as in the paper:
# [dc, Gamma, M_eff, Vb, mass_ratio, u_sdof, d_norm..., d_scaled..., f_push...]
row = [dc, Gamma, M_eff, Vb, mass_ratio, u_sdof] + list(d_norm) + list(d_scaled) + list(f_push)
np.savetxt(os.path.join(outdir, f"pushover_results_{tag}.txt"), np.array([row], float), fmt="%.3f")

save_pushover_layout(
    os.path.join(outdir, f"pushover_layout_{tag}.json"),
    x_d,
    x_stiff_out,
    n_ortho_springs_user,
    n_mains_left,
    n_mains_right,
    x_center,
)

np.savez(
    os.path.join(outdir, f"{tag}.npz"),
    nT=nT,
    nL=nL,
    dc=dc,
    Gamma=Gamma,
    M_eff=M_eff,
    u_sdof=u_sdof,
    nDOF=len(d0),
    nOrth=len(x_ortho_user),
    d_norm=np.array(d_norm),
    stiff_mask=np.array(stiff_mask)
)