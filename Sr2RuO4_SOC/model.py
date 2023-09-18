import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, iOmega_n, inverse, MeshBrZone, MeshProduct
from triqs.lattice import BravaisLattice, BrillouinZone
from triqs.operators import c, c_dag, n
from triqs.operators.util import h_int_kanamori, U_matrix_kanamori
from itertools import product
from numpy import matrix, array, diag, pi
import numpy.linalg as linalg

from tight_binding_model import *

# ==== System Parameters ====
beta = 25.                      # Inverse temperature
mu = 5.3938                     # Chemical potential

U = 2.3                         # Density-density interaction
J = 0.4                         # Hunds coupling
SOC = 0.1                       # Spin-orbit coupling

n_iw = int(10 * beta)           # The number of positive Matsubara frequencies
n_k = 16                        # The number of k-points per dimension

block_names = ['up', 'dn']       # The spins
orb_names = [0, 1, 2]           # The orbitals
idx_lst = list(range(len(block_names) * len(orb_names)))
gf_struct = [('bl', idx_lst)]

TBL = tight_binding_model(lambda_soc=SOC)   # The Tight-Binding Lattice
TBL.bz = BrillouinZone(TBL.bl)
n_idx = len(idx_lst)


# ==== Local Hamiltonian ====
c_dag_vec = matrix([[c_dag('bl', idx) for idx in idx_lst]])
c_vec =     matrix([[c('bl', idx)] for idx in idx_lst])

h_0_mat = TBL._hop[(0,0,0)]
h_0 = (c_dag_vec * h_0_mat * c_vec)[0,0]

Umat, Upmat = U_matrix_kanamori(len(orb_names), U_int=U, J_hund=J)
op_map = { (s,o): ('bl',i) for i, (s,o) in enumerate(product(block_names, orb_names)) }
h_int = h_int_kanamori(block_names, len(orb_names), Umat, Upmat, J, off_diag=True, map_operator_structure=op_map)
h_imp = h_0 + h_int


# ==== Non-Interacting Impurity Green function  ====
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
k_mesh = MeshBrZone(TBL.bz, n_k)
k_iw_mesh = MeshProduct(k_mesh, iw_mesh)

G0_k_iw = BlockGf(mesh=k_iw_mesh, gf_struct=gf_struct)
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)

iw_vec = array([iw.value * np.eye(n_idx) for iw in iw_mesh])
k_vec = array([k.value for k in k_mesh])
e_k_vec = TBL.hopping(k_vec.T.copy() / 2. / pi).transpose(2, 0, 1)
mu_mat = mu * np.eye(n_idx)

G0_k_iw['bl'].data[:] = linalg.inv(iw_vec[None,...] + mu_mat[None,None,...] - e_k_vec[::,None,...])
G0_iw['bl'].data[:] = np.sum(G0_k_iw['bl'].data[:], axis=0) / len(k_mesh)


# ==== Hybridization Function ====
Delta = G0_iw.copy()
Delta['bl'] << iOmega_n + mu_mat - h_0_mat - inverse(G0_iw['bl'])
