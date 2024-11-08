import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from itertools import product
from numpy import sign

# ==== System Parameters ====
beta = 5.           # Inverse temperature
mu = 2.             # Chemical potential
U = 5.              # On-site density-density interaction
h = 0.2             # Local magnetic field
Gamma = 1.          # Hybridization energy

block_names = ['up', 'dn']
n_orb = 1

# ==== Local Hamiltonian ====
h_0 = - mu*( n('up',0) + n('dn',0) ) - h*( n('up',0) - n('dn',0) )
h_int = U * n('up',0) * n('dn',0)
h_imp = h_0 + h_int

# ==== Green function structure ====
gf_struct = [ (s, n_orb) for s in block_names ]

# ==== Hybridization Function ====
n_iw = int(10 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
Delta_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
Delta_iw << 0.;

for iw in Delta_iw['up'].mesh:
    Delta_iw['up'][iw] = -1j * Gamma * sign(iw.imag)
    Delta_iw['dn'][iw] = -1j * Gamma * sign(iw.imag)

# ==== Non-Interacting Impurity Green function  ====
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
G0_iw['up'] << inverse(iOmega_n + mu + h - Delta_iw['up'])
G0_iw['dn'] << inverse(iOmega_n + mu - h - Delta_iw['dn'])
