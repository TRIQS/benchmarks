import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import BlockGf, MeshImFreq, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from triqs.operators.util import h_int_kanamori, U_matrix_kanamori
from itertools import product
from numpy import matrix, array, block, diag, eye
from numpy.linalg import inv

# ==== System Parameters ====
beta = 5.                       # Inverse temperature
mu = 0.0                        # Chemical potential
eps = array([0.0, 0.1])         # Impurity site energies
t = 0.0                         # Hopping between impurity sites

eps_bath = array([0.27, -0.4])  # Bath site energies
t_bath = 0.0                    # Hopping between bath sites

U = 1.                          # Density-density interaction
J = 0.0                         # Hunds coupling

block_names = ['up', 'dn']
n_orb = len(eps)
n_orb_bath = len(eps_bath)

# Non-interacting impurity Hamiltonian in matrix representation
h_0_mat = diag(eps - mu) - matrix([[0, t],
                                   [t, 0]])

# Bath Hamiltonian in matrix representation
h_bath_mat = diag(eps_bath) - matrix([[0, t_bath],
                                      [t_bath, 0]])

# Coupling matrix
V_mat = matrix([[1., 1.],
                [1., 1.]])

# ==== Local Hamiltonian ====
c_dag_vec = { s: matrix([[c_dag(s,o) for o in range(n_orb)]]) for s in block_names }
c_vec =     { s: matrix([[c(s,o)] for o in range(n_orb)]) for s in block_names }

h_0 = sum(c_dag_vec[s] * h_0_mat * c_vec[s] for s in block_names)[0,0]

Umat, Upmat = U_matrix_kanamori(n_orb, U_int=U, J_hund=J)
h_int = h_int_kanamori(block_names, n_orb, Umat, Upmat, J, off_diag=True)

h_imp = h_0 + h_int

# ==== Bath & Coupling Hamiltonian ====
c_dag_bath_vec = { s: matrix([[c_dag(s, o) for o in range(n_orb, n_orb + n_orb_bath)]]) for s in block_names }
c_bath_vec =     { s: matrix([[c(s, o)] for o in range(n_orb, n_orb + n_orb_bath)]) for s in block_names }

h_bath = sum(c_dag_bath_vec[s] * h_bath_mat * c_bath_vec[s] for s in block_names)[0,0]
h_coup = sum(c_dag_vec[s] * V_mat * c_bath_vec[s] + c_dag_bath_vec[s] * V_mat.H * c_vec[s] for s in block_names)[0,0]

# ==== Total impurity Hamiltonian ====
h_tot = h_imp + h_coup + h_bath

# ==== Green function structure ====
gf_struct = [ (s, n_orb) for s in block_names ]

# ==== Non-Interacting Impurity Green function  ====
n_iw = int(10 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
h_tot_mat = block([[h_0_mat, V_mat     ],
                   [V_mat.H, h_bath_mat]])
for bl, iw in product(block_names, iw_mesh):
    G0_iw[bl][iw] = inv(iw.value * eye(2*n_orb) - h_tot_mat)[:n_orb, :n_orb]

# ==== Hybridization Function ====
Delta_iw = G0_iw.copy()
Delta_iw['up'] << iOmega_n - h_0_mat - inverse(G0_iw['up'])
Delta_iw['dn'] << iOmega_n - h_0_mat - inverse(G0_iw['dn'])
