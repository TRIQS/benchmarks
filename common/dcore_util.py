#!/usr/bin/env python
"""Conversion utilities between TRIQS and DCore formats."""

import numpy
from dcore.tools import make_block_gf
from dcore._dispatcher import GfImFreq
from triqs.operators.util.extractors import extract_U_dict4

# DCore block name mapping
_bname_to_dcore = {'up': 'up', 'dn': 'down', 'bl': 'ud'}
_bname_from_dcore = {v: k for k, v in _bname_to_dcore.items()}


def _block_size(entry):
    """Block size from gf_struct entry (name, int) or legacy (name, list)."""
    return entry[1] if isinstance(entry[1], int) else len(entry[1])


def convert_to_dcore_format(gf_struct, h_int, G0_iw, beta, n_iw):
    """Convert TRIQS model data to DCore format (dcorelib types, DCore block names)."""
    use_soc = len(gf_struct) == 1
    norb = _block_size(gf_struct[0]) // 2 if use_soc else _block_size(gf_struct[0])

    # Build DCore gf_struct with list indices
    gf_struct_dcore = {_bname_to_dcore[bl]: list(range(_block_size((bl, n))))
                       for bl, n in gf_struct}

    # Copy G0_iw data into dcorelib Green function
    G0_iw_dcore = make_block_gf(GfImFreq, gf_struct_dcore, beta, n_iw)
    for bl in gf_struct:
        G0_iw_dcore[_bname_to_dcore[bl[0]]].data[:] = G0_iw[bl[0]].data[:]

    # Build Coulomb tensor
    if use_soc:
        idx_tr = {('bl', i): i for i in range(2 * norb)}
    else:
        idx_tr = {(b, i): i + ispin * norb
                  for ispin, b in enumerate(['up', 'dn']) for i in range(norb)}
    U_dict = extract_U_dict4(h_int)
    u_mat = numpy.zeros((2 * norb,) * 4, dtype=complex)
    for idx4, v in U_dict.items():
        idx4_ = [idx_tr[idx] for idx in idx4]
        u_mat[idx4_[0], idx4_[1], idx4_[2], idx4_[3]] += v

    return gf_struct_dcore, u_mat, G0_iw_dcore


def convert_to_triqs_BlockGf(G_dcore, gf_struct, beta, n_iw):
    """Convert dcorelib BlockGf (DCore block names) to TRIQS BlockGf."""
    from triqs.gf import Gf, MeshImFreq, BlockGf
    mesh = MeshImFreq(beta, 'Fermion', n_iw)
    bl_list, g_list = [], []
    for bname_dcore, g in G_dcore:
        bl = _bname_from_dcore[bname_dcore]
        g_triqs = Gf(mesh=mesh, target_shape=g.data.shape[1:])
        g_triqs.data[:] = g.data[:]
        bl_list.append(bl)
        g_list.append(g_triqs)
    return BlockGf(name_list=bl_list, block_list=g_list)
