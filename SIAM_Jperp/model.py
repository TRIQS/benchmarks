import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshDLRImFreq, MeshReFreq, MeshReFreqPts, MeshReFreqLog, BlockGf, Omega, iOmega_n, inverse, SemiCircular
from triqs.gf.descriptors import Function
from triqs.operators import c, c_dag, n
from numpy import matrix

# ==== System Parameters ====
beta = 10.0
U = 1.0
mu = U / 4.0
w0 = 1.0
J = 2.0

block_names = ['up', 'dn']
n_orb = 1
broadening = 1e-3

# ==== Operator vectors (needed for ctseg) ====
c_dag_vec = { s: matrix([[c_dag(s,o) for o in range(n_orb)]]) for s in block_names }
c_vec     = { s: matrix([[c(s,o)] for o in range(n_orb)]) for s in block_names }

# ==== Local Hamiltonian ====
h_int = U * n('dn', 0) * n('up', 0)

# ==== Green function structure ====
gf_struct = [(bl, n_orb) for bl in block_names]

# ==== Frequency Meshes ====
n_iw = int(10 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
dlr_wmax = 10.0
dlr_eps = 1e-10
dlr_iw_mesh = MeshDLRImFreq(beta, 'Fermion', dlr_wmax, dlr_eps, True)

# ==== Non-Interacting Impurity Green function and Hybridization ====
def make_g0_and_delta(mesh):
    if type(mesh) in [MeshReFreq, MeshReFreqPts, MeshReFreqLog]:
        z = Omega + 1j * broadening
    else:
        z = iOmega_n
    G0 = BlockGf(mesh=mesh, gf_struct=gf_struct)
    for bl, g_bl in G0:
        g_bl << inverse(z + mu - SemiCircular(1.0))
    Delta = G0.copy()
    for bl, d_bl in Delta:
        d_bl << z + mu - inverse(G0[bl])
    return G0, Delta

G0_iw, Delta_iw = make_g0_and_delta(iw_mesh)
G0_dlr_iw, Delta_dlr = make_g0_and_delta(dlr_iw_mesh)

# ==== Dynamic spin-spin interaction Jperp ====
# Jperp(iw) = 0.5 * J^2 * (1/(w-w0) - 1/(w+w0))
Jperp_func = Function(lambda w: 0.5 * J**2 * (1/(w - w0) - 1/(w + w0)))
