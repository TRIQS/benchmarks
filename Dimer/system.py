execfile('../common/util.py')

from pytriqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from pytriqs.operators import c, c_dag, n
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from itertools import product
from numpy import matrix, array

# ==== System Parameters ====
beta = 5.           # Inverse temperature
mu = 0.25           # Chemical potential
U = 1.              # Density-density interaction for opposite spins
Up = 0.3            # Density-density interaction for equal spins
J = 0.5             # Hunds coupling
t = 1.              # Hopping between impurity sites
epsilon = matrix([[0.2,0.1],[0.1,0.2]]) # Bath state energy

# ==== Local Hamiltonian ====
h_0 = - mu*( n('up',0) + n('dn',0) + n('up',1) + n('dn',1) ) \
      + t* ( c_dag('up',0) * c('up',1) + c_dag('up',1) * c('up',0) + \
             c_dag('dn',0) * c('dn',1) + c_dag('dn',1) * c('dn',0) )

h_int = h_int_kanamori(['up','dn'],[0,1],
                        array([[0,Up-3*J],[Up-3*J,0]]), # Interaction for equal spins
                        array([[U,U-2*J],[U-2*J,U]]),   # Interaction for opposite spins
                        J,True)

h_loc = h_0 + h_int

# ==== Bath & Coupling hamiltonian ====
h_bath, h_coup = 0, 0
for sig in ['up','dn']:
    h_coup += c_dag(sig,0) * c(sig,2) + c_dag(sig,2) * c(sig,0)
    h_coup += c_dag(sig,1) * c(sig,3) + c_dag(sig,3) * c(sig,1)
    h_bath += epsilon[0,0] * n(sig,2) + epsilon[1,1] * n(sig,3)
    h_bath += epsilon[0,1] * c_dag(sig,2) * c(sig,3) + epsilon[1,0] * c_dag(sig,3) * c(sig,2)

# ==== Total impurity hamiltonian and fundamental operators ====
h_imp = h_loc + h_coup + h_bath
fundamental_operators = [ c(spin,i) for spin, i in product(['up','dn'],range(4)) ]

# ==== Green function structure ====
gf_struct = [ ['up',[0,1]], ['dn',[0,1]] ]

# ==== Hybridization Function ====
n_iw = 10
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
Delta = Gf_from_struct(mesh=iw_mesh, struct=gf_struct)
Delta << inverse(iOmega_n - epsilon)

# ==== Non-Interacting Impurity Green function  ====
G0_iw = Gf_from_struct(mesh=iw_mesh, struct=gf_struct)
G0_iw['up'] << inverse(iOmega_n + mu - t * matrix([[0,1],[1,0]]) - Delta['up']) # FIXME Should work for BlockGf
G0_iw['dn'] << inverse(iOmega_n + mu - t * matrix([[0,1],[1,0]]) - Delta['dn'])
