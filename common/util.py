import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from triqs.gf import *
from triqs.operators import c, c_dag

def get_fundamental_operators(op):
    """Extract unique annihilation operators from a many-body operator."""
    idx_lst = []
    for term, val in op:
        for has_dagger, (bl, orb) in term:
            if (bl, orb) not in idx_lst:
                idx_lst.append((bl, orb))
    return [c(bl, orb) for bl, orb in idx_lst]


def Block2Gf_from_struct(mesh, struct):
    """Create a Block2Gf container with spin-block structure.

    For chi3: mesh = MeshProduct(iw_mesh, iw_mesh), target_shape = (bl1_size, bl1_size, bl2_size, bl2_size)
    For chi2: mesh = iW_mesh, target_shape = (bl1_size, bl1_size, bl2_size, bl2_size)
    """
    bl_lst = [bl for bl, bl_size in struct]
    G_lst = []
    for bl1, bl1_size in struct:
        lst = []
        for bl2, bl2_size in struct:
            lst.append(Gf(mesh=mesh, target_shape=(bl1_size, bl1_size, bl2_size, bl2_size)))
        G_lst.append(lst)
    return Block2Gf(name_list1=bl_lst, name_list2=bl_lst, block_list=G_lst)


def Block2Gf_from_fourier2D(G_tau, n_iw, struct):
    """2D Fourier transform of a Block2Gf from tau x tau to iw x iw."""
    bl1, G1 = next(iter(G_tau))
    tau_meshes = [G1.mesh[0], G1.mesh[1]]
    beta = tau_meshes[0].beta
    iw_meshes = [MeshImFreq(beta, tau_mesh.statistic, n_iw) for tau_mesh in tau_meshes]

    G_iw = Block2Gf_from_struct(mesh=MeshProduct(*iw_meshes), struct=struct)
    temp = Block2Gf_from_struct(mesh=MeshProduct(iw_meshes[0], tau_meshes[1]), struct=struct)

    for bl, G_tau_bl in G_tau:
        for tau in tau_meshes[1]:
            temp[bl][:, tau] << Fourier(G_tau_bl[:, tau])
        for iw in iw_meshes[0]:
            G_iw[bl][iw, :] << Fourier(temp[bl][iw, :])

    return G_iw


def kronecker(iw, iwp):
    """Kronecker delta for Matsubara frequencies."""
    return iw == iwp


def make_gf_dlr_iw(*G_iw_list, dlr_eps=1e-10, wmax_init=1.0, wmax_max=200.0, fit_eps_factor=0.1):
    """Find optimal DLR wmax and return DLR representations of the input Green functions.

    Takes one or more BlockGf on MeshImFreq (typically G0_iw, Delta_iw).
    Returns (*G_dlr_iw_list, dlr_wmax) -- DLR BlockGfs + the determined wmax.

    The wmax is determined by gradually increasing it until the DLR round-trip
    of the first argument (G0_iw) reproduces the original within dlr_eps.
    A tighter eps (dlr_eps * fit_eps_factor) is used for DLR basis construction.
    """
    import numpy as np
    from triqs.gf import Gf, MeshDLRImFreq, make_gf_dlr, make_gf_imfreq

    G0_iw = G_iw_list[0]
    mesh = G0_iw.mesh
    beta = mesh.beta
    statistic = str(mesh.statistic)
    n_iw = len(mesh) // 2
    fit_eps = dlr_eps * fit_eps_factor

    wmax = wmax_init
    while wmax <= wmax_max:
        dlr_mesh = MeshDLRImFreq(beta, statistic, wmax, fit_eps, True)

        # Check DLR round-trip accuracy on G0_iw
        max_err = 0.0
        for bl, g_bl in G0_iw:
            g_dlr_iw = Gf(mesh=dlr_mesh, target_shape=g_bl.target_shape)
            for w in dlr_mesh:
                g_dlr_iw[w] = g_bl(w)
            g_dlr = make_gf_dlr(g_dlr_iw)
            g_rec = make_gf_imfreq(g_dlr, n_iw)
            err = np.max(np.abs(g_rec.data - g_bl.data))
            max_err = max(max_err, err)

        if max_err < dlr_eps:
            # Build DLR GFs for all inputs on the converged mesh
            results = []
            for G_iw in G_iw_list:
                G_dlr = BlockGf(mesh=dlr_mesh, gf_struct=[(bl, g.target_shape[0]) for bl, g in G_iw])
                for bl, g_bl in G_iw:
                    for w in dlr_mesh:
                        G_dlr[bl][w] = g_bl(w)
                results.append(G_dlr)
            results.append(wmax)
            return tuple(results)

        wmax *= 1.5

    raise RuntimeError(f"make_gf_dlr_iw: could not find suitable wmax <= {wmax_max} for dlr_eps={dlr_eps}")
