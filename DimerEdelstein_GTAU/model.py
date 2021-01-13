import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from triqs.operators.util import h_int_kanamori, U_matrix_kanamori
from itertools import product
from numpy import matrix, array, block, diag, eye
from numpy.linalg import inv

from math import sqrt

# ==== System Parameters ====
beta = 32.                      # Inverse temperature
mu = 0.0                        # Chemical potential
eps = array([0.0, 0.0])         # Impurity site energies
t = 0.0                         # Hopping between impurity sites

eps_bath = array([-2.3, 2.3])  # Bath site energies
t_bath = 0.0                       # Hopping between bath sites

U = 2                         # Density-density interaction
J = 0.2                       # Hunds coupling

spin_names = ['up', 'dn']
orb_names  = [0, 1]
n_orb = len(orb_names)

# Non-interacting impurity hamiltonian in matrix representation
h_0_mat = diag(eps - mu) - matrix([[0, t],
                                   [t, 0]])

# Bath hamiltonian in matrix representation
h_bath_mat = diag(eps_bath) - matrix([[0, t_bath],
                                      [t_bath, 0]])

# Coupling matrix
r = 0.5
a = 100. # 20.0
V_mat = sqrt(a) * 1. / 2. * matrix([[sqrt(1+r) + sqrt(1-r), sqrt(1+r) - sqrt(1-r)],
                          [sqrt(1+r) - sqrt(1-r), sqrt(1+r) + sqrt(1-r)]])
# V_mat * V_mat.T = a * [[ 1,  r ],
#                        [ r,  1 ]]

# ==== Local Hamiltonian ====
c_dag_vec = { s: matrix([[c_dag(s,o) for o in orb_names]]) for s in spin_names }
c_vec =     { s: matrix([[c(s,o)] for o in orb_names]) for s in spin_names }

h_0 = sum(c_dag_vec[s] * h_0_mat * c_vec[s] for s in spin_names)[0,0]

## h_int_kanamori includes also pair-hopping terms! Not considered in Edelstein Paper
#Umat, Upmat = U_matrix_kanamori(n_orb, U_int=U, J_hund=J)
#h_int = h_int_kanamori(spin_names, orb_names, Umat, Upmat, J, off_diag=True)
h_int = U * (n('up', 0) * n('dn', 0) + n('up', 1) * n('dn', 1)) \
     + (U - 2*J) * (n('up', 0) * n('dn', 1) + n('up', 1) * n('dn', 0)) \
     + (U - 3*J) * (n('up', 0) * n('up', 1) + n('dn', 0) * n('dn', 1)) \
     + J * c_dag('up', 0) * c_dag('dn', 1) * c('dn', 0) * c('up', 1) \
     + J * c_dag('up', 1) * c_dag('dn', 0) * c('dn', 1) * c('up', 0)
#h_int = 0.01 * h_int

h_imp = h_0 + h_int

# ==== Bath & Coupling hamiltonian ====
orb_bath_names = [n_orb + o for o in orb_names]
c_dag_bath_vec = { s: matrix([[c_dag(s, o) for o in orb_bath_names]]) for s in spin_names }
c_bath_vec =     { s: matrix([[c(s, o)] for o in orb_bath_names]) for s in spin_names }

h_bath = sum(c_dag_bath_vec[s] * h_bath_mat * c_bath_vec[s] for s in spin_names)[0,0]
h_coup = sum(c_dag_vec[s] * V_mat * c_bath_vec[s] + c_dag_bath_vec[s] * V_mat * c_vec[s] for s in spin_names)[0,0] # FIXME Adjoint

# ==== Total impurity hamiltonian ====
h_tot = h_imp + h_coup + h_bath

# ==== Green function structure ====
# gf_struct = [ [s, orb_names] for s in spin_names ]
gf_struct = [ [s, n_orb] for s in spin_names ]

# ==== Non-Interacting Impurity Green function  ====
n_iw = int(10 * beta)
#n_tau = 2 * n_iw + 1
n_tau = 20
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
tau_mesh = MeshImTime(beta, 'Fermion', n_tau)
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
h_tot_mat = block([[h_0_mat, V_mat     ],
                   [V_mat.H, h_bath_mat]])
for bl, iw in product(spin_names, iw_mesh):
    G0_iw[bl][iw] = inv(iw.value * eye(2*n_orb) - h_tot_mat)[:n_orb, :n_orb]

# ==== Hybridization Function ====
Delta = G0_iw.copy()
Delta['up'] << iOmega_n - h_0_mat - inverse(G0_iw['up'])
Delta['dn'] << iOmega_n - h_0_mat - inverse(G0_iw['dn'])
