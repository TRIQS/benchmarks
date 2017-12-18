import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from triqs.gf import *
from triqs.operators import c, c_dag, n, dagger

# Get a list of all annihilation operators from a many-body operators
def get_fundamental_operators(op):
    idx_lst = []
    for term, val in op:
        for has_dagger, (bl, orb) in term:
            if not idx_lst.count([bl, orb]):
                idx_lst.append([bl,orb])
    return [c(bl, orb) for bl, orb in idx_lst]

# print on master node
def mpi_print(arg):
    if mpi.is_master_node():
        print(arg)

def Block2Gf_from_struct(mesh, struct):
    
    ## Without block structure
    #if not isinstance(struct[0], list):
    #    return Gf(mesh=mesh, indices=struct)

    bl_lst = [bl for bl, bl_size in struct]

    # With block structure
    G_lst = []
    for bl1, bl1_size in struct:
        lst = []
        for bl2, bl2_size in struct:
            lst.append(Gf(mesh=mesh, target_shape=(bl1_size, bl1_size, bl2_size, bl2_size)))
        G_lst.append(lst)

    return Block2Gf(name_list1=bl_lst, name_list2=bl_lst, block_list=G_lst)

# Not Working, tail problems...
def Block2Gf_from_fourier2D(G_tau, n_iw, struct):
    
    bl1, G1 = next(iter(G_tau))

    tau_meshes = [G1.mesh[0], G1.mesh[1]]
    beta = tau_meshes[0].beta
    iw_meshes = [MeshImFreq(beta, tau_mesh.statistic, n_iw) for tau_mesh in tau_meshes]

    #FIXME struct should be member of Block2Gf
    G_iw = Block2Gf_from_struct(mesh=MeshProduct(*iw_meshes), struct=struct)
    temp = Block2Gf_from_struct(mesh=MeshProduct(iw_meshes[0], tau_meshes[1]), struct=struct)

    for bl, G_tau_bl in G_tau:
        for tau in tau_meshes[1]:
            temp[bl][:,tau] << Fourier(G_tau_bl[:,tau])
        for iw in iw_meshes[0]:
            G_iw[bl][iw,:] << Fourier(temp[bl][iw,:])

    return G_iw
    
def kronecker(iw, iwp):
    return iw == iwp

def norm_inf(G):
    import numpy

    if isinstance(G, (BlockGf, Block2Gf)):
        return max(norm_inf(G_bl) for bl, G_bl in G)

    elif isinstance(G, Gf):
        return norm_inf(G.data)

    return numpy.amax(numpy.abs(G))
