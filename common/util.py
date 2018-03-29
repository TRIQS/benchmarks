from pytriqs.gf import *
from pytriqs.operators import c, c_dag, n, dagger

# Build a block Green function for a given gf_struct
def BlockGf_from_struct(mesh, struct):
    # Without block structure
    if not isinstance(struct[0], list):
        return Gf(mesh=mesh, indices=struct)

    # With block structure
    G_lst = []
    for bl, idx_lst in struct:
        G_lst.append(Gf(mesh=mesh, indices=idx_lst))
    return BlockGf(name_list=[bl[0] for bl in struct], block_list=G_lst)

# Get a list of all annihilation operators from a many-body operators
def get_fundamental_operators(op):
    idx_lst = []
    for term, val in op:
        for has_dagger, (bl, orb) in term:
            if not idx_lst.count([bl, orb]):
                idx_lst.append([bl,orb])
    return [c(bl, orb) for bl, orb in idx_lst]
