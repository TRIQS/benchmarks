from pytriqs.gf import *
from pytriqs.operators import c, c_dag, n, dagger

# Get a list of all annihilation operators from a many-body operators
def get_fundamental_operators(op):
    idx_lst = []
    for term, val in op:
        for has_dagger, (bl, orb) in term:
            if not idx_lst.count([bl, orb]):
                idx_lst.append([bl,orb])
    return [c(bl, orb) for bl, orb in idx_lst]
