import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from triqs.gf import *
from triqs.gf.tools import fit_legendre
from triqs.operators import c, c_dag, n, dagger

# Get a list of all annihilation operators from a many-body operators
def get_fundamental_operators(op):
    idx_lst = []
    for term, val in op:
        for has_dagger, (bl, orb) in term:
            if not idx_lst.count([bl, orb]):
                idx_lst.append([bl,orb])
    return [c(bl, orb) for bl, orb in idx_lst]

def fit_G_l(Gt):
    G_l = fit_legendre(Gt, 5)

    for n_l in range(7, 41, 2):
        G_l_new = fit_legendre(Gt, n_l)
        if abs(G_l_new['up'][0,0].data[n_l-3]) < abs(G_l_new['up'][0,0].data[n_l-1]):
            print("Optimal number of legendre coefficients " + str(n_l - 2))
            return G_l
        else:
            G_l = G_l_new

    return G_l
