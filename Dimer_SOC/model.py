import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from pytriqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from pytriqs.operators import c, c_dag, n
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from itertools import product
from numpy import matrix, array, diag, eye
from numpy.linalg import inv

# ==== System Parameters ====
beta = 5.                       # Inverse temperature
mu = 0.25                       # Chemical potential
eps = array([0.0, 0.1, 0.0, 0.1])   # Impurity site energies
t = 1.                          # Hopping between impurity sites
a = 1j                          # Spin-orbit coupling

eps_bath = array([0.2, 0.15, 0.2, 0.15])   # Bath site energies
t_bath = 0.1                    # Hopping between bath sites

U = 1.                          # Density-density interaction for opposite spins
Up = 0.3                        # Density-density interaction for equal spins

up_0, up_1, dn_0, dn_1 = 0, 1, 2, 3
orb_names  = [up_0, up_1, dn_0, dn_1]

# Non-interacting impurity hamiltonian in matrix representation
h_0_mat = diag(eps - mu) - matrix([[0,    t+a,  0,    a  ],
                                   [t-a,  0,   -a,    0  ],
                                   [0,    a,    0,    t+a],
                                   [-a,    0,    t-a,  0  ]])

# Bath hamiltonian in matrix representation
h_bath_mat = diag(eps_bath) - matrix([[0,       t_bath, 0,      0     ],
                                      [t_bath,  0,      0,      0     ],
                                      [0,       0,      0,      t_bath],
                                      [0,       0,      t_bath, 0     ]])

# Coupling matrix
V_mat = matrix([[1., 0., 0, 0],
                [0., 1., 0, 0],
                [0., 0., 1, 0],
                [0., 0., 0, 1]])

# ==== Local Hamiltonian ====
c_dag_vec = matrix([[c_dag('bl',o) for o in orb_names]])
c_vec = matrix([[c('bl',o)] for o in orb_names])

h_0 = (c_dag_vec * h_0_mat * c_vec)[0,0]

h_int = U * n('bl', up_0) * n('bl', dn_0) + \
        U * n('bl', up_1) * n('bl', dn_1) + \
        U * n('bl', up_0) * n('bl', dn_1) + \
        U * n('bl', up_1) * n('bl', dn_0) + \
        Up * n('bl', up_0) * n('bl', up_1) + \
        Up * n('bl', dn_0) * n('bl', dn_1)

h_loc = h_0 + h_int

# ==== Bath & Coupling hamiltonian ====
orb_bath_names = ['b_' + str(o) for o in orb_names]
c_dag_bath_vec = matrix([[c_dag('bl', o) for o in orb_bath_names]])
c_bath_vec =     matrix([[c('bl', o)] for o in orb_bath_names])

h_bath = (c_dag_bath_vec * h_bath_mat * c_bath_vec)[0,0]
h_coup = (c_dag_vec * V_mat * c_bath_vec + c_dag_bath_vec * V_mat * c_vec)[0,0] # FIXME Adjoint

# ==== Total impurity hamiltonian ====
h_imp = h_loc + h_coup + h_bath

# ==== Green function structure ====
gf_struct = [ ['bl', orb_names] ]

# ==== Hybridization Function ====
n_iw = int(10 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
Delta = BlockGf(mesh=iw_mesh, gf_struct=gf_struct)
# FIXME Delta['bl'] << V_mat * inverse(iOmega_n - h_bath_mat) * V_mat.transpose()
for iw in iw_mesh:
    Delta['bl'][iw] = V_mat * inv(iw.value * eye(len(orb_names)) - h_bath_mat) * V_mat.transpose()

# ==== Non-Interacting Impurity Green function  ====
G0_iw = Delta.copy()
G0_iw['bl'] << inverse(iOmega_n - h_0_mat - Delta['bl'])
