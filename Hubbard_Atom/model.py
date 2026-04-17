r"""# Hubbard Atom

A single atomic level with on-site Coulomb repulsion $U$, chemical potential $\mu$,
and local magnetic (Zeeman) field $h$. No bath — the impurity is isolated, so the
Green function of the atom can be written down analytically.

$$
H = -\mu\,(n_\uparrow + n_\downarrow)
    - h\,(n_\uparrow - n_\downarrow)
    + U\, n_\uparrow n_\downarrow
$$
"""
import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshReFreq, MeshReFreqPts, MeshReFreqLog, Omega, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from itertools import product

# ==== System Parameters ====
beta = 5.           # Inverse temperature
mu = 2.             # Chemical potential
U = 5.              # On-site density-density interaction
h = 0.2             # Local magnetic field

block_names = ['up', 'dn']
n_orb = 1
n_orb_bath = 0

# ==== Local Hamiltonian ====
h_0 = - mu*( n('up',0) + n('dn',0) ) - h*( n('up',0) - n('dn',0) )
h_int = U * n('up',0) * n('dn',0)
h_imp = h_0 + h_int

# ==== Bath & Coupling Hamiltonian ====
h_bath, h_coup = 0, 0

# ==== Total impurity hamiltonian and fundamental operators ====
h_tot = h_imp + h_coup + h_bath

# ==== Green function structure ====
gf_struct = [ (s, n_orb) for s in block_names ]

# ==== Frequency Meshes ====
n_iw = int(10 * beta)
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)

# ==== Real-Frequency Mesh (for models with discrete spectra) ====
n_w = 3001
w_window = (-10, 10)
w_mesh = MeshReFreq(window=w_window, n_w=n_w)
broadening = 1e-3

# ==== Non-Interacting Impurity Green function and target Hybridization ====
def make_g0_and_delta(mesh):
    if type(mesh) in [MeshReFreq, MeshReFreqPts, MeshReFreqLog]:
        z = Omega + 1j * broadening
    else:
        z = iOmega_n
    Delta = BlockGf(mesh=mesh, gf_struct=gf_struct)
    Delta << 0.0
    G0 = BlockGf(mesh=mesh, gf_struct=gf_struct)
    G0['up'] << inverse(z + mu + h - Delta['up'])
    G0['dn'] << inverse(z + mu - h - Delta['dn'])
    return G0, Delta

G0_iw, Delta_iw = make_g0_and_delta(iw_mesh)

# ==== DLR Green functions (auto-determined wmax) ====
dlr_eps = 1e-10
G0_dlr_iw, Delta_dlr, dlr_wmax = make_gf_dlr_iw(G0_iw, Delta_iw, dlr_eps=dlr_eps)
dlr_iw_mesh = G0_dlr_iw.mesh
