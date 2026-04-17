r"""# Trimer — Kanamori + discrete bath

Three-orbital (three-site) impurity with inter-site hopping $t$, a Kanamori
interaction $(U, U', J)$, and a three-level discrete bath — a direct
generalisation of `Dimer`. The non-interacting part couples the three sites
via the hopping matrix, and each site hybridises with its own bath level.

See `Dimer` for the Kanamori form of $H_{\mathrm{int}}$.
"""
import sys, os
sys.path.append(os.getcwd() + '/../common')
from util import *

from triqs.gf import Gf, MeshImFreq, MeshReFreq, MeshReFreqPts, MeshReFreqLog, Omega, iOmega_n, inverse
from triqs.operators import c, c_dag, n
from triqs.operators.util import h_int_kanamori, U_matrix_kanamori
from itertools import product
from numpy import matrix, array, block, diag, eye
from numpy.linalg import inv

# ==== System Parameters ====
beta = 5.                               # Inverse temperature
mu = 0.0                                # Chemical potential
eps = array([0.0, 0.1, 0.2])            # Impurity site energies
t = 0.2                                 # Hopping between impurity sites

eps_bath = array([0.27, -0.4, 0.15])    # Bath site energies
t_bath = 0.0                            # Hopping between bath sites

U = 1.                                  # Density-density interaction
J = 0.2                                 # Hunds coupling

block_names = ['up', 'dn']
n_orb = len(eps)
n_orb_bath = len(eps_bath)

# Non-interacting impurity hamiltonian in matrix representation
h_0_mat = diag(eps - mu) - matrix([[0, t, t],
                                   [t, 0, t],
                                   [t, t, 0]])

# Bath hamiltonian in matrix representation
h_bath_mat = diag(eps_bath) - matrix([[0,       t_bath, t_bath  ],
                                      [t_bath,  0,      t_bath  ],
                                      [t_bath,  t_bath, 0       ]])

# Coupling matrix
V_mat = matrix([[1., 1., 1.],
                [1., 1., 1.],
                [1., 1., 1.]])

# ==== Local Hamiltonian ====
c_dag_vec = { s: matrix([[c_dag(s,o) for o in range(n_orb)]]) for s in block_names }
c_vec =     { s: matrix([[c(s,o)] for o in range(n_orb)]) for s in block_names }

h_0 = sum(c_dag_vec[s] * h_0_mat * c_vec[s] for s in block_names)[0,0]

Umat, Upmat = U_matrix_kanamori(n_orb, U_int=U, J_hund=J)
h_int = h_int_kanamori(block_names, n_orb, Umat, Upmat, J, off_diag=True)

h_imp = h_0 + h_int

# ==== Bath & Coupling hamiltonian ====
c_dag_bath_vec = { s: matrix([[c_dag(s, o) for o in range(n_orb, n_orb + n_orb_bath)]]) for s in block_names }
c_bath_vec =     { s: matrix([[c(s, o)] for o in range(n_orb, n_orb + n_orb_bath)]) for s in block_names }

h_bath = sum(c_dag_bath_vec[s] * h_bath_mat * c_bath_vec[s] for s in block_names)[0,0]
h_coup = sum(c_dag_vec[s] * V_mat * c_bath_vec[s] + c_dag_bath_vec[s] * V_mat * c_vec[s] for s in block_names)[0,0] # FIXME Adjoint

# ==== Total impurity hamiltonian ====
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

# ==== Non-Interacting Impurity Green function and Hybridization ====
h_tot_mat = block([[h_0_mat, V_mat     ],
                   [V_mat.H, h_bath_mat]])

def make_g0_and_delta(mesh):
    if type(mesh) in [MeshReFreq, MeshReFreqPts, MeshReFreqLog]:
        z = Omega + 1j * broadening
    else:
        z = iOmega_n
    real_freq = type(mesh) in [MeshReFreq, MeshReFreqPts, MeshReFreqLog]
    G0 = BlockGf(mesh=mesh, gf_struct=gf_struct)
    for bl, w in product(block_names, mesh):
        zv = w.value + 1j * broadening if real_freq else w.value
        G0[bl][w] = inv(zv * eye(2*n_orb) - h_tot_mat)[:n_orb, :n_orb]
    Delta = G0.copy()
    Delta['up'] << z - h_0_mat - inverse(G0['up'])
    Delta['dn'] << z - h_0_mat - inverse(G0['dn'])
    return G0, Delta

G0_iw, Delta_iw = make_g0_and_delta(iw_mesh)

# ==== DLR Green functions (auto-determined wmax) ====
dlr_eps = 1e-10
G0_dlr_iw, Delta_dlr, dlr_wmax = make_gf_dlr_iw(G0_iw, Delta_iw, dlr_eps=dlr_eps)
dlr_iw_mesh = G0_dlr_iw.mesh
