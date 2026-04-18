import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshDLRImFreq, BlockGf, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from numpy import array, matrix, eye, zeros

# ==== System Parameters ====
beta = 80.0                     # Inverse temperature
U = 2.0                         # Hubbard interaction
t = 1.0                         # Hopping
mu = U / 2                      # Chemical potential (half-filling)

Nx, Ny = 4, 4                   # Cluster dimensions
n_orb = Nx * Ny                 # 16 sites in 4x4 cluster
n_orb_bath = 0
n_iw = int(5 * beta)            # Matsubara frequencies

block_names = ['up', 'dn']

# ==== Hopping matrix (4x4 cluster with PBC) ====
hopping = zeros((Nx, Ny, Nx, Ny))
for d in range(Nx):
    hopping[d, :, (d + 1) % Nx, :] = t * eye(Ny)
    hopping[(d + 1) % Nx, :, d, :] = t * eye(Ny)
for d in range(Ny):
    hopping[:, d, :, (d + 1) % Ny] = t * eye(Nx)
    hopping[:, (d + 1) % Ny, :, d] = t * eye(Nx)
hopping = hopping.reshape(n_orb, n_orb)

h_0_mat = -mu * eye(n_orb) - hopping

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
dlr_wmax = 8.0
dlr_eps = 1e-8
dlr_iw_mesh = MeshDLRImFreq(beta, 'Fermion', dlr_wmax, dlr_eps, True)

# ==== Non-Interacting Green function and Hybridization (zero for isolated cluster) ====
def make_gf(mesh):
    G0 = BlockGf(mesh=mesh, gf_struct=gf_struct)
    for bl, g_bl in G0:
        g_bl << inverse(iOmega_n - h_0_mat)
    Delta = BlockGf(mesh=mesh, gf_struct=gf_struct)
    Delta << 0.0
    return G0, Delta

G0_iw, Delta_iw = make_gf(iw_mesh)
G0_dlr_iw, Delta_dlr = make_gf(dlr_iw_mesh)
