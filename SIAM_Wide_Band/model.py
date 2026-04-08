import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshDLRImFreq, MeshReFreq, MeshReFreqPts, MeshReFreqLog, Omega, iOmega_n, inverse
from triqs.gf.descriptors import Function
from triqs.operators import c, c_dag, n
from itertools import product
from numpy import sign, matrix

# ==== System Parameters ====
beta = 5.           # Inverse temperature
mu = 2.             # Chemical potential
U = 5.              # On-site density-density interaction
h = 0.2             # Local magnetic field
Gamma = 1.          # Hybridization energy

block_names = ['up', 'dn']
n_orb = 1

# ==== Operator vectors (needed for ctseg) ====
c_dag_vec = { s: matrix([[c_dag(s,o) for o in range(n_orb)]]) for s in block_names }
c_vec     = { s: matrix([[c(s,o)] for o in range(n_orb)]) for s in block_names }

# ==== Local Hamiltonian ====
h_0 = - mu*( n('up',0) + n('dn',0) ) - h*( n('up',0) - n('dn',0) )
h_int = U * n('up',0) * n('dn',0)
h_imp = h_0 + h_int

# ==== Green function structure ====
gf_struct = [ (s, n_orb) for s in block_names ]

# ==== Frequency Meshes ====
n_iw = int(10 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
dlr_wmax = 2 * U
dlr_eps = 1e-10
dlr_iw_mesh = MeshDLRImFreq(beta, 'Fermion', dlr_wmax, dlr_eps, True)

# ==== Non-Interacting Impurity Green function and Hybridization ====
def make_g0_and_delta(mesh):
    Delta = BlockGf(mesh=mesh, gf_struct=gf_struct)
    if type(mesh) in [MeshReFreq, MeshReFreqPts, MeshReFreqLog]:
        z = Omega
        Delta << -1j * Gamma
    else:
        z = iOmega_n
        Delta << Function(lambda w: -1j * Gamma * sign(w.value.imag))
    G0 = BlockGf(mesh=mesh, gf_struct=gf_struct)
    G0['up'] << inverse(z + mu + h - Delta['up'])
    G0['dn'] << inverse(z + mu - h - Delta['dn'])
    return G0, Delta

G0_iw, Delta_iw = make_g0_and_delta(iw_mesh)
G0_dlr_iw, Delta_dlr = make_g0_and_delta(dlr_iw_mesh)
