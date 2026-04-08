import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshReFreq, MeshReFreqPts, MeshReFreqLog, Omega, iOmega_n, inverse, MeshBrZone, MeshProduct
from triqs.lattice import BravaisLattice, BrillouinZone
from triqs.operators import c, c_dag, n
from triqs.operators.util import h_int_kanamori, U_matrix_kanamori
from itertools import product
import numpy as np
from numpy import matrix, array, diag, pi
import numpy.linalg as linalg

from triqs.lattice.utils import TB_from_wannier90

# order is xy_up, xz_up, yz_up, xy_dn, xz_dn, yz_dn
def lambda_matrix(lam_xy, lam_z):
    lam_loc = np.zeros((6,6), dtype=complex)
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
n_idx = len(idx_lst)
gf_struct = [('bl', n_idx)]

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
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
k_mesh = MeshBrZone(TBL.bz, n_k)

e_k_vec = TBL.fourier(k_mesh)
mu_mat = mu * np.eye(n_idx)
broadening = 1e-3

# ==== Matsubara Green function and Hybridization ====
def make_g0_and_delta(mesh):
    real_freq = type(mesh) in [MeshReFreq, MeshReFreqPts, MeshReFreqLog]
    if real_freq:
        z = Omega + 1j * broadening
        z_vec = array([(w.value + 1j * broadening) * np.eye(n_idx) for w in mesh])
    else:
        z = iOmega_n
        z_vec = array([w.value * np.eye(n_idx) for w in mesh])
    G0 = BlockGf(mesh=mesh, gf_struct=gf_struct)
    G0_k = linalg.inv(z_vec[None,...] + mu_mat[None,None,...] - e_k_vec.data[::,None,...])
    G0['bl'].data[:] = np.sum(G0_k, axis=0) / len(k_mesh)
    Delta = G0.copy()
    Delta['bl'] << z + mu_mat - h_0_mat - inverse(G0['bl'])
    return G0, Delta

G0_iw, Delta_iw = make_g0_and_delta(iw_mesh)

# ==== DLR Green functions (auto-determined wmax) ====
dlr_eps = 1e-10
G0_dlr_iw, Delta_dlr, dlr_wmax = make_gf_dlr_iw(G0_iw, Delta_iw, dlr_eps=dlr_eps)
dlr_iw_mesh = G0_dlr_iw.mesh

# ==== k-resolved Green function ====
k_iw_mesh = MeshProduct(k_mesh, iw_mesh)
G0_k_iw = BlockGf(mesh=k_iw_mesh, gf_struct=gf_struct)
iw_vec = array([iw.value * np.eye(n_idx) for iw in iw_mesh])
G0_k_iw['bl'].data[:] = linalg.inv(iw_vec[None,...] + mu_mat[None,None,...] - e_k_vec.data[::,None,...])
