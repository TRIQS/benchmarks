import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, iOmega_n, inverse, MeshBrZone, MeshProduct
from triqs.lattice import BravaisLattice, BrillouinZone
from triqs.operators import c, c_dag, n
from triqs.operators.util import h_int_kanamori, U_matrix_kanamori
from itertools import product
#import numpy as np
from numpy import matrix, array, diag, pi
import numpy.linalg as linalg

from triqs.lattice.utils import TB_from_wannier90

# ==== System Parameters ====
beta = 25.                     # Inverse temperature
mu = 5.3938                    # Chemical potential

U = 2.3                         # Density-density interaction
J = 0.4                         # Hunds coupling

n_iw = int(10 * beta)           # The number of positive Matsubara frequencies
n_k = 16                        # The number of k-points per dimension

block_names = ['up', 'dn']       # The spins
n_orb = 3

paths = [os.getcwd(), os.path.dirname(__file__)]
for p in paths:
    if os.path.isfile(p + '/w2w_hr.dat'):
        path = p
        break
TBL = TB_from_wannier90(seed='/w2w', path=path)
TBL.bz = BrillouinZone(TBL.bl)


# ==== Local Hamiltonian ====
c_dag_vec = { s: matrix([[c_dag(s,o) for o in range(n_orb)]]) for s in block_names }
c_vec =     { s: matrix([[c(s,o)] for o in range(n_orb)]) for s in block_names }

h_0_mat = TBL.hoppings[(0,0,0)][0:n_orb,0:n_orb]
h_0 = sum(c_dag_vec[s] * h_0_mat * c_vec[s] for s in block_names)[0,0]

Umat, Upmat = U_matrix_kanamori(n_orb, U_int=U, J_hund=J)
h_int = h_int_kanamori(block_names, n_orb, Umat, Upmat, J, off_diag=True)

h_imp = h_0 + h_int


# ==== Non-Interacting Impurity Green function  ====
gf_struct = [ (s,n_orb) for s in block_names ]

iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
k_mesh = MeshBrZone(TBL.bz, n_k)
k_iw_mesh = MeshProduct(k_mesh, iw_mesh)

G0_k_iw = BlockGf(mesh=k_iw_mesh, gf_struct=gf_struct)
G0_iw = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)

iw_vec = array([iw.value * np.eye(n_orb) for iw in iw_mesh])
e_k = TBL.fourier(k_mesh)[0:n_orb,0:n_orb]
mu_mat = mu * np.eye(n_orb)

for s in block_names:
    G0_k_iw[s].data[:] = linalg.inv(iw_vec[None,...] + mu_mat[None,None,...] - e_k.data[::,None,...])
    G0_iw[s].data[:] = np.sum(G0_k_iw[s].data[:], axis=0) / len(k_mesh)


# ==== Hybridization Function ====
Delta = G0_iw.copy()
Delta['up'] << iOmega_n + mu_mat - h_0_mat - inverse(G0_iw['up'])
Delta['dn'] << iOmega_n + mu_mat - h_0_mat - inverse(G0_iw['dn'])
