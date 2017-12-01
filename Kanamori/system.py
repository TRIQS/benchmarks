execfile('../common/util.py')

from pytriqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from pytriqs.operators import c, c_dag, n
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from itertools import product
import numpy as np

# ==== System Parameters ====
beta = 5.           # Inverse temperature
mu = 2.             # Chemical potential
U = 5.              # On-site density-density interaction
J = 0.2             # Hunds coupling
V = 1.0             # Bath coupling
epsilon = 2.3       # Bath state energy

# ==== Local Hamiltonian ====
J = 0.

h_0 = - mu*( n('up',0) + n('dn',0) + n('up',1) + n('dn',1) )
h_int = h_int_kanamori(['up','dn'],[0,1],
                        # np.array([[0,U-3*J],[U-3*J,0]]),
                        np.array([[0,0],[0,0]]),
                        # np.array([[U,U-2*J],[U-2*J,U]]),
                        np.array([[U,0],[0,U]]),
                        J,True)
h_loc = h_0 + h_int

# ==== Green function structure ====
gf_struct = [ ['up',[0,1]], ['dn',[0,1]] ]

# ==== Hybridization Function ====
n_iw = 20
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
Delta = Gf_from_struct(mesh=iw_mesh, struct=gf_struct)
Delta << (V**2) * inverse(iOmega_n - epsilon) + (V**2) * inverse(iOmega_n + epsilon);

# ==== Non-Interacting Impurity Green function  ====
G0_iw = Gf_from_struct(mesh=iw_mesh, struct=gf_struct)
G0_iw['up'] << 1.0 * inverse(iOmega_n + mu - Delta['up']) # FIXME Should work for BlockGf
G0_iw['dn'] << 1.0 * inverse(iOmega_n + mu - Delta['dn'])
