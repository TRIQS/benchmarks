import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from itertools import product

# ==== System Parameters ====
# Parameters from arXiv:2405.06716 (SIAM with discrete bath)
beta = 160.         # Inverse temperature
U = 5.              # On-site density-density interaction
mu = U / 2.         # Chemical potential (half-filling)
h = 0.0             # Local magnetic field
E = [0.0, -U, U]    # Bath-site energies (bandwidth 2U)
V = [1.0, 1.0, 1.0] # Couplings to Bath-sites

block_names = ['up', 'dn']
n_orb = 1
n_orb_bath = len(E)

# ==== Local Hamiltonian ====
h_0 = - mu*( n('up',0) + n('dn',0) ) - h*( n('up',0) - n('dn',0) )
h_int = U * n('up',0) * n('dn',0)
h_imp = h_0 + h_int

# ==== Bath & Coupling Hamiltonian ====
h_bath, h_coup = 0, 0
for i, (E_i, V_i) in enumerate(zip(E, V)):
    for sig in ['up','dn']:
        h_bath += E_i * n(sig, n_orb + i)
        h_coup += V_i * (c_dag(sig,0) * c(sig, n_orb + i) + c_dag(sig, n_orb + i) * c(sig,0))

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

# ==== Hybridization Functions ====
Delta = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
Delta << sum([V_i*V_i * inverse(iOmega_n - E_i) for V_i, E_i in zip(V, E)]);
Delta_dlr = BlockGf(mesh=dlr_iw_mesh, gf_struct=gf_struct)
Delta_dlr << sum([V_i*V_i * inverse(iOmega_n - E_i) for V_i, E_i in zip(V, E)]);

# ==== Non-Interacting Impurity Green functions ====
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
G0_iw['up'] << inverse(iOmega_n + mu + h - Delta['up'])
G0_iw['dn'] << inverse(iOmega_n + mu - h - Delta['dn'])
G0_dlr_iw = BlockGf(mesh=dlr_iw_mesh, gf_struct=gf_struct)
G0_dlr_iw['up'] << inverse(iOmega_n + mu + h - Delta_dlr['up'])
G0_dlr_iw['dn'] << inverse(iOmega_n + mu - h - Delta_dlr['dn'])
