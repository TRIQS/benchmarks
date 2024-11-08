import sys
import os

import numpy as np

sys.path.append(os.getcwd() + '/../common')

from triqs.gf import BlockGf, MeshImFreq, MeshReFreq, Omega, iOmega_n, inverse, MeshBrZone, MeshProduct
from triqs.lattice import BrillouinZone
from triqs.operators import c, c_dag
from triqs.operators.util import h_int_kanamori, U_matrix_kanamori
import numpy.linalg as linalg

from triqs.lattice.utils import TB_from_wannier90

# ==== System Parameters ====
beta = 25.0  # Inverse temperature
mu = 12.79  # Chemical potential
eta = 2e-1

U = 3.6  # Density-density interaction

n_iw = int(10 * beta)  # The number of positive Matsubara frequencies
n_w = 3001
window = [-10, 10]
n_k = 41  # The number of k-points per dimensio

block_names = ['up', 'dn']  # The spins
n_orb = 1

paths = [os.getcwd(), os.path.dirname(__file__)]
for p in paths:
    if os.path.isfile(p + '/lco_hr.dat'):
        path = p
        break
TBL = TB_from_wannier90(seed='/lco', path=path)
TBL.bz = BrillouinZone(TBL.bl)


# ==== Local Hamiltonian ====
c_dag_vec = {s: np.matrix([[c_dag(s, o) for o in range(n_orb)]]) for s in block_names}
c_vec = {s: np.matrix([[c(s, o)] for o in range(n_orb)]) for s in block_names}

h_0_mat = TBL.hoppings[(0, 0, 0)][0:n_orb, 0:n_orb]
h_0 = sum(c_dag_vec[s] * h_0_mat * c_vec[s] for s in block_names)[0, 0]

Umat, Upmat = U_matrix_kanamori(n_orb, U_int=U, J_hund=0.0)
h_int = h_int_kanamori(block_names, n_orb, Umat, Upmat, 0.0, off_diag=True)

h_imp = h_0 + h_int

# ==== Non-Interacting Impurity Green function  ====
gf_struct = [(s, n_orb) for s in block_names]

iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
w_mesh = MeshReFreq(window, n_w)
k_mesh = MeshBrZone(TBL.bz, n_k)
k_iw_mesh = MeshProduct(k_mesh, iw_mesh)
k_w_mesh = MeshProduct(k_mesh, w_mesh)

G0_k_iw = BlockGf(mesh=k_iw_mesh, gf_struct=gf_struct)
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
G0_k_w = BlockGf(mesh=k_w_mesh, gf_struct=gf_struct)
G0_w = BlockGf(mesh=w_mesh, gf_struct=gf_struct)

iw_vec = np.array([iw.value * np.eye(n_orb) for iw in iw_mesh])
w_vec = np.array([w.value * np.eye(n_orb) for w in w_mesh])
e_k = TBL.fourier(k_mesh)[0:n_orb, 0:n_orb]
mu_mat = mu * np.eye(n_orb)
dc_mat = 1.9 * np.eye(n_orb)
eta_mat = 1j * eta * np.eye(n_orb)

for s in block_names:
    G0_k_iw[s].data[:] = linalg.inv(iw_vec[None, ...] + mu_mat[None, None, ...] - e_k.data[::, None, ...])
    G0_iw[s].data[:] = np.sum(G0_k_iw[s].data[:], axis=0) / len(k_mesh)

    G0_k_w[s].data[:] = linalg.inv(w_vec[None, ...] + mu_mat[None, None, ...] - e_k.data[::, None, ...] + eta_mat[None, None, ...])
    G0_w[s].data[:] = np.sum(G0_k_w[s].data[:], axis=0) / len(k_mesh)


# ==== Hybridization Function ====
Delta_iw = G0_iw.copy()
Delta_iw['up'] << iOmega_n + mu_mat - h_0_mat - inverse(G0_iw['up'])
Delta_iw['dn'] << iOmega_n + mu_mat - h_0_mat - inverse(G0_iw['dn'])

# ReFreq Delta
Delta_w = G0_w.copy()
Delta_w['up'] << Omega + mu_mat + eta_mat - h_0_mat - inverse(G0_w['up'])
Delta_w['dn'] << Omega + mu_mat + eta_mat - h_0_mat - inverse(G0_w['dn'])
