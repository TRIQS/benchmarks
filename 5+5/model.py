import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from pytriqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from pytriqs.operators import c, c_dag, n
from pytriqs.operators.util.U_matrix import cubic_names, U_matrix
from pytriqs.operators.util.hamiltonians import h_int_slater
from itertools import product
from numpy import matrix, array

# ==== System Parameters ====
beta = 5.           # Inverse temperature
mu = 26.            # Chemical potential
U = 2.              # Density-density interaction for opposite spins
J = 0.5             # Hunds coupling
F0 = U
F2 = J*(14.0/(1.0 + 0.63))
F4 = F2*0.63
U_mat = U_matrix(2,[F0,F2,F4],basis='cubic')
orb_names = list(cubic_names(2))

# Hybridization function parameters
delta_params={"xy"      : {'V':0.2,'e':-0.2},
              "yz"      : {'V':0.2,'e':-0.15},
              "z^2"     : {'V':0.2,'e':-0.1},
              "xz"      : {'V':0.2,'e':0.05},
              "x^2-y^2" : {'V':0.2,'e':0.4}}

atomic_levels = {('up', 'xy')       : -0.2,
                 ('dn', 'xy')       : -0.2,
                 ('up', 'yz')       : -0.15,
                 ('dn', 'yz')       : -0.15,
                 ('up', 'z^2')      : -0.1,
                 ('dn', 'z^2')      : -0.1,
                 ('up', 'xz')       : 0.05,
                 ('dn', 'xz')       : 0.05,
                 ('up', 'x^2-y^2')  : 0.4,
                 ('dn', 'x^2-y^2')  : 0.4}

# ==== Local Hamiltonian ====
# h_0 = ... 

h_int = h_int_slater(['up','dn'], orb_names, U_mat, True)

# h_loc = h_0 + h_int

# ==== Bath & Coupling hamiltonian ====
# h_bath, h_coup = 0, 0
# h_bath = ... 
# h_coup = ... 

# ==== Total impurity hamiltonian and fundamental operators ====
# h_imp = h_loc + h_coup + h_bath
# fundamental_operators = [ c(spin,i) for spin, i in product(['up','dn'],range(4)) ]

# ==== Green function structure ====
gf_struct = [['up', orb_names ], ['dn', orb_names ] ]

# ==== Hybridization Function & Non-Interacting Impurity Green function====
n_iw = 10
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
Delta = Gf_from_struct(mesh=iw_mesh, struct=gf_struct)
G0_iw = Gf_from_struct(mesh=iw_mesh, struct=gf_struct)

G0_iw << inverse(iOmega_n)

for spin, orbital in product(['up','dn'],orb_names):
    V = delta_params[orbital]['V']
    e = delta_params[orbital]['e']

    Delta[spin][orbital,orbital] << (V**2) * inverse(iOmega_n - e)
    G0_iw[spin][orbital,orbital] << inverse(iOmega_n +mu - atomic_levels[(spin,orbital)] - Delta[spin][orbital,orbital])
