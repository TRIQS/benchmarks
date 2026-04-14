import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshDLRImFreq, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from itertools import product

# ==== System Parameters ====
# Parameters from arXiv:2405.06716 (half-filled Hubbard atom)
beta = 160.         # Inverse temperature
U = 1.              # On-site density-density interaction
mu = U / 2.         # Chemical potential (half-filling condition)
h = 0.              # No magnetic field

block_names = ['up', 'dn']
n_orb = 1
n_orb_bath = 0

# ==== Local Hamiltonian ====
h_0 = - mu*( n('up',0) + n('dn',0) ) - h*( n('up',0) - n('dn',0) )
h_int = U * n('up',0) * n('dn',0)
h_imp = h_0 + h_int

# ==== Bath & Coupling Hamiltonian ====
h_bath, h_coup = 0, 0

# ==== Total impurity hamiltonian and fundamental operators ====
h_tot = h_imp + h_coup + h_bath

# ==== Green function structure ====
gf_struct = [ (s, n_orb) for s in block_names ]

# ==== Frequency Meshes ====
n_iw = int(5 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
dlr_wmax = 2*U
dlr_eps = 1e-8
dlr_iw_mesh = MeshDLRImFreq(beta, 'Fermion', dlr_wmax, dlr_eps, True)

# ==== Non-Interacting Impurity Green function and target Hybridization ====
def make_gf(mesh):
    Delta = BlockGf(mesh=mesh, gf_struct=gf_struct)
    Delta << 0.0
    G0 = BlockGf(mesh=mesh, gf_struct=gf_struct)
    G0['up'] << inverse(iOmega_n + mu + h - Delta['up'])
    G0['dn'] << inverse(iOmega_n + mu - h - Delta['dn'])
    return G0, Delta

G0_iw, Delta_iw = make_gf(iw_mesh)
G0_dlr_iw, Delta_dlr = make_gf(dlr_iw_mesh)
