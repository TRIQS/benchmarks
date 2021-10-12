import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from itertools import product

from numpy import matrix

# ==== System Parameters ====
beta = 5.           # Inverse temperature
U = 2.              # On-site density-density interaction
mu = U / 2          # Chemical potential
h = 0.0             # Local magnetic field
d = 0.0             # Local pairing field

block_names = ['bl']
n_orb = 2
n_orb_bath = 0

# ==== Local Hamiltonian in Nambu Basis Psi = (d_up, ddag_dn) ====
# n_up <-> n_0
# n_dn <-> (1 - n_1)
h0_mat = matrix([[-mu-h+U,         d    ],
                 [ d.conjugate(),  mu-h ]])
h_0 = (matrix([[c_dag('bl',0), c_dag('bl',1)]]) * h0_mat * matrix([[c('bl',0)], [c('bl',1)]]))[0,0]

h_int = -U * n('bl',0) * n('bl',1)
h_imp = h_0 + h_int

# ==== Bath & Coupling Hamiltonian ====
h_bath, h_coup = 0, 0

# ==== Total impurity hamiltonian and fundamental operators ====
h_tot = h_imp + h_coup + h_bath

# ==== Green function structure ====
gf_struct = [ (s, n_orb) for s in block_names ]

# ==== Hybridization Function ====
n_iw = int(10 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
Delta = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
Delta << 0.0;

# ==== Non-Interacting Impurity Green function  ====
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
G0_iw['bl'] << inverse(iOmega_n - h0_mat - Delta['bl'])

print("h_0: ", h_0)
print("h_0: ", h_0 + 1)
print("G0_iw[0,0].density(): ", G0_iw['bl'][0,0].density())
print("G0_iw[1,1].density(): ", G0_iw['bl'][1,1].density())
