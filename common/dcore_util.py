#!/usr/bin/env python
from __future__ import print_function

from dcore.tools import make_block_gf, raise_if_mpi_imported

import sys, os, time
import numpy
from itertools import product

import pytriqs
from pytriqs.gf import GfImFreq
from pytriqs.operators.util.extractors import extract_U_dict4

raise_if_mpi_imported()

def convert_to_dcore_format(gf_struct, h_int, G0_iw, beta, n_iw):
     """
     Convert input data to DCore format:
         * gf_struct is a dict.
         * Allowed block names are "ud", "up", "down".
     """
     # --------- Convert gf_struct ----------
     num_blocks = len(gf_struct)
     for b in gf_struct:
         assert b[0] in ['up', 'dn', 'bl']
     bname_tr = {'up':'up', 'dn':'down', 'bl':'ud'}
     gf_struct_dcore = {bname_tr[b[0]] : b[1] for b in gf_struct}
     gf_struct_dcore_list = [(bname_tr[b[0]], b[1]) for b in gf_struct]
     norb = len(gf_struct[0][1])//2 if num_blocks == 1 else len(gf_struct[0][1])

     if num_blocks == 1:
         idx_tr = {('bl', i): i for i in range(2*norb)}
         G0_iw_dcore = make_block_gf(GfImFreq, gf_struct_dcore, beta, n_iw)
         G0_iw_dcore['ud'] = G0_iw['bl']
     else:
         idx_tr = {(b, i): i+ispin*norb for ispin, b in enumerate(['up', 'dn']) for i in range(norb)}
         G0_iw_dcore = make_block_gf(GfImFreq, gf_struct_dcore, beta, n_iw)
         G0_iw_dcore['up'] = G0_iw['up']
         G0_iw_dcore['down'] = G0_iw['dn']

     # --------- Construct Coulomb tensor from h_int ----------
     U_dict = extract_U_dict4(h_int)
     u_mat = numpy.zeros((2*norb,)*4, dtype=complex)
     for idx4, v in U_dict.items():
         print("U_dict: ", idx4, v)
         idx4_ = map(lambda idx: idx_tr[idx], idx4)
         u_mat[idx4_[0], idx4_[1], idx4_[2], idx4_[3]] += v
         print(idx4, idx4_, u_mat[idx4_[0], idx4_[1], idx4_[2], idx4_[3]])

     return norb, gf_struct_dcore, u_mat, G0_iw_dcore

def convert_to_triqs_bname(G, gf_struct, beta, n_iw):
    gf_struct_dict = {b[0] : b[1] for b in gf_struct}
    G_copy = make_block_gf(GfImFreq, gf_struct_dict, beta, n_iw)
    bname_trans = {'up':'up', 'down':'dn', 'ud':'bl'}
    for bname, g in G:
        G_copy[bname_trans[bname]] = g
    return G_copy
