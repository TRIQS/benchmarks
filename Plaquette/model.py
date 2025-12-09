import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, BlockGf, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from numpy import array, matrix

# ==== System Parameters ====
beta = 25.0                     # Inverse temperature
U = 2.0                         # Hubbard interaction
t = 1.0                         # Hopping
mu = U / 2                      # Chemical potential (half-filling)

n_orb = 4                       # 4 sites in 2x2 plaquette
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

# ==== Green function structure ====
gf_struct = [(s, n_orb) for s in block_names]

# ==== Non-Interacting Green function ====
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
for bl, g_bl in G0_iw:
    g_bl << inverse(iOmega_n - h_0_mat)

# ==== Hybridization Function (zero for isolated cluster) ====
Delta_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
Delta_iw << 0.0
