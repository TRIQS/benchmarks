import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from pytriqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from pytriqs.operators import c, c_dag, n
from itertools import product
from numpy import sign

# ==== System Parameters ====
beta = 5.           # Inverse temperature
mu = 2.             # Chemical potential
U = 5.              # On-site density-density interaction
h = 0.2             # Local magnetic field
Gamma = 1.          # Hybridization energy

# ==== Local Hamiltonian ====
h_0 = - mu*( n('up',0) + n('dn',0) ) - h*( n('up',0) - n('dn',0) )
h_int = U * n('up',0) * n('dn',0)
h_loc = h_0 + h_int

# ==== Green function structure ====
gf_struct = [ ['up',[0]], ['dn',[0]] ]

# ==== Hybridization Function ====
n_iw = 20
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
Delta = BlockGf_from_struct(mesh=iw_mesh, struct=gf_struct)
Delta << 0.; 

for iw in Delta['up'].mesh:
    Delta['up'][iw] = -1j * Gamma * sign(iw.imag)
    Delta['dn'][iw] = -1j * Gamma * sign(iw.imag)

# ==== Non-Interacting Impurity Green function  ====
G0_iw = BlockGf_from_struct(mesh=iw_mesh, struct=gf_struct)
G0_iw['up'] << inverse(iOmega_n + mu + h) # FIXME Set tails explicitly
G0_iw['dn'] << inverse(iOmega_n + mu - h) # FIXME Set tails explicitly

for iw in iw_mesh:
    G0_iw['up'][iw] = 1.0 / (iw + mu + h - Delta['up'][iw])
    G0_iw['dn'][iw] = 1.0 / (iw + mu - h - Delta['dn'][iw])
