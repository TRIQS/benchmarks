import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, iOmega_n, inverse, MeshBrZone, MeshProduct
from triqs.lattice import BravaisLattice, BrillouinZone
from triqs.operators import c, c_dag, n
from triqs.operators.util import h_int_kanamori, U_matrix_kanamori
from itertools import product
import numpy as np
from numpy import matrix, array, diag, pi
import numpy.linalg as linalg

from triqs.lattice.utils import TB_from_wannier90

# order is xy_up, xz_up, yz_up, xy_dn, xz_dn, yz_dn 
# could be replaced by triqs_dft_tools/converters/wannier90.py: generate_local_so_matrix_t2g()
# but different orbital order (xz_up, xz_dn, yz_up, yz_dn, xy_up, xy_dn)
def lambda_matrix(lam_xy, lam_z):
    lam_loc = np.zeros((6,6),dtype=complex)
    lam_loc[0,4] =  1j*lam_xy/2.0
    lam_loc[0,5] =     lam_xy/2.0
    lam_loc[1,2] =  1j*lam_z/2.0
    lam_loc[1,3] = -1j*lam_xy/2.0
    lam_loc[2,3] =    -lam_xy/2.0
    lam_loc[4,5] = -1j*lam_z/2.0
    lam_loc = lam_loc + np.transpose(np.conjugate(lam_loc))
    return lam_loc

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
n_orb = len(idx_lst)

paths = [os.getcwd(), os.path.dirname(__file__)]
for p in paths:
    if os.path.isfile(p + '/w2w_hr.dat'):
        path = p
        break
TBL = TB_from_wannier90(seed='/w2w', path=path, extend_to_spin=True, add_local=lambda_matrix(SOC, SOC))
TBL.bz = BrillouinZone(TBL.bl)


# ==== Local Hamiltonian ====
c_dag_vec = matrix([[c_dag('bl', idx) for idx in idx_lst]])
c_vec =     matrix([[c('bl', idx)] for idx in idx_lst])

h_0_mat = TBL.hoppings[(0,0,0)]
h_0 = (c_dag_vec * h_0_mat * c_vec)[0,0]

Umat, Upmat = U_matrix_kanamori(len(orb_names), U_int=U, J_hund=J)
op_map = { (s,o): ('bl',i) for i, (s,o) in enumerate(product(block_names, orb_names)) }
h_int = h_int_kanamori(block_names, len(orb_names), Umat, Upmat, J, off_diag=True, map_operator_structure=op_map)

h_imp = h_0 + h_int


# ==== Non-Interacting Impurity Green function  ====
gf_struct = [('bl', len(idx_lst))]

iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
k_mesh = MeshBrZone(TBL.bz, n_k)
k_iw_mesh = MeshProduct(k_mesh, iw_mesh)

G0_k_iw = BlockGf(mesh=k_iw_mesh, gf_struct=gf_struct)
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)

iw_vec = array([iw.value * np.eye(n_orb) for iw in iw_mesh])
e_k = TBL.fourier(k_mesh)[0:n_orb,0:n_orb]
mu_mat = mu * np.eye(n_orb)

G0_k_iw['bl'].data[:] = linalg.inv(iw_vec[None,...] + mu_mat[None,None,...] - e_k.data[::,None,...])
G0_iw['bl'].data[:] = np.sum(G0_k_iw['bl'].data[:], axis=0) / len(k_mesh)


# ==== Hybridization Function ====
Delta_iw = G0_iw.copy()
Delta_iw['bl'] << iOmega_n + mu_mat - h_0_mat - inverse(G0_iw['bl'])
