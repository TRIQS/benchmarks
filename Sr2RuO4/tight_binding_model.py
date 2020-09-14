
# ----------------------------------------------------------------------

import os
import copy
import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.lattice.tight_binding import TBLattice, energies_on_bz_path

# ----------------------------------------------------------------------

from triqs_tprf.wannier90 import parse_hopping_from_wannier90_hr_dat
from triqs_tprf.wannier90 import parse_lattice_vectors_from_wannier90_wout
from triqs_tprf.wannier90 import parse_band_structure_from_wannier90_band_dat

# ----------------------------------------------------------------------
def extend_wannier90_to_spin(hopping, num_wann):

    hopping_spin = {}
    for key, value in list(hopping.items()):
        #hopping_spin[key] = np.kron(value, np.eye(2)) # spin fastest idx
        hopping_spin[key] = np.kron(np.eye(2), value) # orbital fastest idx

    return hopping_spin, 2 * num_wann

# ----------------------------------------------------------------------
def tight_binding_model(crystal_field=0., lambda_soc=0.): 

    paths = [os.getcwd(), os.path.dirname(__file__)]
    for p in paths:
        if os.path.isfile(p + '/w2w_hr.dat'):
            path = p
            break

    # -- Read Wannier90 results

    hopping, num_wann = parse_hopping_from_wannier90_hr_dat(path + '/w2w_hr.dat')
    orbital_names = [str(i) for i in range(num_wann)]
    units = parse_lattice_vectors_from_wannier90_wout(path + '/w2w.wout')

    # -- Extend to spinful model from non-spin polarized Wannier90 result
    hopping_spin, num_wann_spin = extend_wannier90_to_spin(hopping, num_wann)
    #orbital_names_spin = ["".join(reversed(tup)) for tup in itertools.product(orbital_names, ['up_', 'do_'])]
    orbital_names_spin = ["".join(tup) for tup in itertools.product(['up_', 'do_'], orbital_names)]

    # ------------------------------------------------------------------
    # order is xy_up, xz_up, yz_up, xy_dn, xz_dn, yz_dn 

    def lambda_matrix_pavarini(lam_xy, lam_z): # according to https://arxiv.org/pdf/1612.03060.pdf
        lam_loc = np.zeros((6,6),dtype=complex)
        lam_loc[0,5] =     lam_xy/2.0
        lam_loc[0,4] =  1j*lam_xy/2.0
        lam_loc[1,2] = -1j*lam_z/2.0
        lam_loc[2,3] =    -lam_xy/2.0
        lam_loc[1,3] = -1j*lam_xy/2.0
        lam_loc[4,5] =  1j*lam_z/2.0
        lam_loc = lam_loc + np.transpose(np.conjugate(lam_loc))
        return lam_loc

    def lambda_matrix(lam_xy, lam_z):
        lam_loc = np.zeros((6,6),dtype=complex)
        lam_loc[0,4] =  1j*lam_xy/2.0
        lam_loc[0,5] =     lam_xy/2.0
        lam_loc[1,2] =  1j*lam_z/2.0
        lam_loc[1,3] = -1j*lam_xy/2.0
        lam_loc[2,3] =    -lam_xy/2.0
        lam_loc[4,5] = -1j*lam_z/2.0
        lam_loc = lam_loc + np.transpose(np.conjugate(lam_loc))
        return lam_loc

    def cf_matrix(cf_xy, cf_z):
        cf_loc = np.zeros((6,6),dtype=complex)
        cf_loc[0,0] =    cf_xy
        cf_loc[1,1] =    cf_z
        cf_loc[2,2] =    cf_z
        cf_loc[3,3] =    cf_xy
        cf_loc[4,4] =    cf_z
        cf_loc[5,5] =    cf_z
        return cf_loc

    #lam_loc = lambda_matrix(0.11,0.11)
    #print 'Lambda Matrix: \n', lam_loc

    #lam_loc_pavarini = lambda_matrix_pavarini(0.11,0.11)
    #print 'Lambda Matrix Pavarini: \n', lam_loc_pavarini

    #d_lam_loc = lambda_matrix(0.20,0.20)
    #print 'Lambda + dLambda Matrix: \n', d_lam_loc

    #cf_loc = cf_matrix(-0.01,+0.01)

    cf_loc = cf_matrix(-crystal_field, +crystal_field)
    d_lam_loc = lambda_matrix(lambda_soc, lambda_soc)    

    # add soc and cf terms
    hopping = copy.deepcopy(hopping_spin)
    hopping[(0,0,0)] += d_lam_loc+ cf_loc

    # ------------------------------------------------------------------

    #print '--> tight binding model'

    tb_lattice = TBLattice(
        units = units,
        hopping = hopping,
        orbital_positions = [(0,0,0)]*num_wann_spin,
        orbital_names = orbital_names_spin,
        )

    tb_lattice.lambda_soc = lambda_soc
    tb_lattice.crystal_field = crystal_field
    
    return tb_lattice
