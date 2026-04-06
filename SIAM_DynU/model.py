import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshDLRImFreq, BlockGf, iOmega_n, inverse, SemiCircular
from triqs.gf.descriptors import Function
from triqs.operators import c, c_dag, n
from numpy import matrix

# ==== System Parameters ====
beta = 10.0
U = 1.0
mu = U / 4.0
w0 = 1.0
D_coupling = 0.5

block_names = ['dn', 'up']
n_orb = 1

# ==== Operator vectors (needed for ctseg) ====
c_dag_vec = { s: matrix([[c_dag(s,o) for o in range(n_orb)]]) for s in block_names }
c_vec     = { s: matrix([[c(s,o)] for o in range(n_orb)]) for s in block_names }

# ==== Local Hamiltonian ====
h_int = U * n('dn', 0) * n('up', 0)

# ==== Green function structure ====
gf_struct = [(bl, n_orb) for bl in block_names]

# ==== Non-Interacting Impurity Green function ====
n_iw = int(10 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
for bl, g_bl in G0_iw:
    g_bl << inverse(iOmega_n + mu - SemiCircular(1.0))

# ==== DLR Green function ====
dlr_wmax = 10.0
dlr_eps = 1e-10
dlr_iw_mesh = MeshDLRImFreq(beta, 'Fermion', dlr_wmax, dlr_eps, True)
G0_dlr_iw = BlockGf(mesh=dlr_iw_mesh, gf_struct=gf_struct)
for bl, g_bl in G0_dlr_iw:
    g_bl << inverse(iOmega_n + mu - SemiCircular(1.0))

# ==== Dynamic density-density interaction D0 ====
# D0(iw) = D_coupling^2 * (1/(w-w0) - 1/(w+w0)) for ('up','dn') and ('dn','up') blocks
D0_func = Function(lambda w: D_coupling**2 * (1/(w - w0) - 1/(w + w0)))
D0_block_pairs = [('up', 'dn'), ('dn', 'up')]
