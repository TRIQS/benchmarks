import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshDLRImFreq, MeshReFreq, MeshReFreqPts, MeshReFreqLog, BlockGf, Omega, iOmega_n, inverse
from triqs.gf.descriptors import Function
from triqs.operators import c, c_dag, n
from numpy import array, sign, eye, matrix

# ==== System Parameters ====
beta = 25.0                     # Inverse temperature
U = 2.0                         # Hubbard interaction
t = 1.0                         # Hopping
mu = U / 2                      # Chemical potential (half-filling)
Gamma = 0.1                     # Hybridization strength (wide-band limit)

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

# ==== Total Hamiltonian (isolated cluster + wide-band bath) ====
h_tot = h_imp

# ==== Green function structure ====
gf_struct = [(s, n_orb) for s in block_names]

# ==== Frequency Meshes ====
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
dlr_wmax = 2 * U
dlr_eps = 1e-10
dlr_iw_mesh = MeshDLRImFreq(beta, 'Fermion', dlr_wmax, dlr_eps, True)

# ==== Non-Interacting Green function and Hybridization (wide-band limit) ====
def make_g0_and_delta(mesh):
    Delta = BlockGf(mesh=mesh, gf_struct=gf_struct)
    if type(mesh) in [MeshReFreq, MeshReFreqPts, MeshReFreqLog]:
        z = Omega
        Delta << -1j * Gamma
    else:
        z = iOmega_n
        Delta << Function(lambda w: -1j * Gamma * sign(w.value.imag))
    G0 = BlockGf(mesh=mesh, gf_struct=gf_struct)
    for bl in block_names:
        G0[bl] << inverse(z - h_0_mat - Delta[bl])
    return G0, Delta

G0_iw, Delta_iw = make_g0_and_delta(iw_mesh)
G0_dlr_iw, Delta_dlr = make_g0_and_delta(dlr_iw_mesh)
