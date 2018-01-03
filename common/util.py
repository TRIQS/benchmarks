from pytriqs.gf import *

def BlockGf_from_struct(mesh, struct):
    # Without block structure
    if not isinstance(struct[0], list):
        return Gf(mesh=mesh, indices=struct)

    # With block structure
    G_lst = []
    for bl, idx_lst in struct:
        G_lst.append(Gf(mesh=mesh, indices=idx_lst))
    return BlockGf(name_list=[bl[0] for bl in struct], block_list=G_lst)
