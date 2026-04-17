r"""# Plaquette — Isolated 2×2 Hubbard cluster

Four-site $2 \times 2$ cluster with periodic boundary conditions, on-site
Hubbard interaction, and **no bath**:

$$
H = -t \sum_{\langle ij\rangle,\sigma}
        (c^\dagger_{i\sigma} c_{j\sigma} + \mathrm{H.c.})
    - \mu \sum_{i\sigma} n_{i\sigma}
    + U \sum_i n_{i\uparrow} n_{i\downarrow},
$$

at half-filling ($\mu = U/2$). Useful as a finite-size cluster benchmark for
cluster DMFT-style solvers.
"""
import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshReFreq, MeshReFreqPts, MeshReFreqLog, BlockGf, Omega, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from numpy import array, matrix

# ==== System Parameters ====
beta = 25.0                     # Inverse temperature
U = 2.0                         # Hubbard interaction
t = 1.0                         # Hopping
mu = U / 2                      # Chemical potential (half-filling)

n_orb = 4                       # 4 sites in 2x2 plaquette
n_orb_bath = 0
n_iw = int(10 * beta)           # Matsubara frequencies

block_names = ['up', 'dn']

# ==== Hopping matrix (2x2 plaquette with PBC) ====
# Site numbering: 0-1
#                 2-3
# Hopping: 0↔1, 1↔3, 3↔2, 2↔0
h_0_mat = -array([
    [mu,  t, t, 0],
    [ t, mu, 0, t],
    [ t,  0, mu, t],
    [ 0,  t, t, mu],
])

# ==== Interaction Hamiltonian ====
h_int = sum(U * n('up', i) * n('dn', i) for i in range(n_orb))

# ==== Local Hamiltonian (quadratic part) ====
c_dag_vec = {s: matrix([[c_dag(s, o) for o in range(n_orb)]]) for s in block_names}
c_vec = {s: matrix([[c(s, o)] for o in range(n_orb)]) for s in block_names}
h_0 = sum(c_dag_vec[s] * h_0_mat * c_vec[s] for s in block_names)[0, 0]

h_imp = h_0 + h_int

# ==== Total Hamiltonian (isolated cluster, no bath) ====
h_tot = h_imp

# ==== Green function structure ====
gf_struct = [(s, n_orb) for s in block_names]

# ==== Frequency Meshes ====
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)

# ==== Real-Frequency Mesh (for models with discrete spectra) ====
n_w = 3001
w_window = (-10, 10)
w_mesh = MeshReFreq(window=w_window, n_w=n_w)
broadening = 1e-3

# ==== Non-Interacting Green function and Hybridization (zero for isolated cluster) ====
def make_g0_and_delta(mesh):
    if type(mesh) in [MeshReFreq, MeshReFreqPts, MeshReFreqLog]:
        z = Omega + 1j * broadening
    else:
        z = iOmega_n
    G0 = BlockGf(mesh=mesh, gf_struct=gf_struct)
    for bl, g_bl in G0:
        g_bl << inverse(z - h_0_mat)
    Delta = BlockGf(mesh=mesh, gf_struct=gf_struct)
    Delta << 0.0
    return G0, Delta

G0_iw, Delta_iw = make_g0_and_delta(iw_mesh)

# ==== DLR Green functions (auto-determined wmax) ====
dlr_eps = 1e-10
G0_dlr_iw, Delta_dlr, dlr_wmax = make_gf_dlr_iw(G0_iw, Delta_iw, dlr_eps=dlr_eps)
dlr_iw_mesh = G0_dlr_iw.mesh
