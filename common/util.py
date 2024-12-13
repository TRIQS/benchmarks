import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from triqs.gf import *
from triqs.operators import c, c_dag, n, dagger

import numpy as np

# Get a list of all annihilation operators from a many-body operators
def get_fundamental_operators(op):
    idx_lst = []
    for term, val in op:
        for has_dagger, (bl, orb) in term:
            if not idx_lst.count([bl, orb]):
                idx_lst.append([bl,orb])
    return [c(bl, orb) for bl, orb in idx_lst]

# Construct a logarithmic mesh; returns a numpy array
def log_mesh(mesh_max, mesh_min, mesh_ratio):
    mesh = np.array([])

    w = mesh_max
    while w > mesh_min:
        mesh = np.append(mesh, w)
        w = w / mesh_ratio

    mesh = np.concatenate((-mesh, np.flip(mesh)))

    return mesh