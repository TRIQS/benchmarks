import sys, os
sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/../../common")
from model import *
import util

from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from pomerol2triqs import PomerolED

# --------- Construct the ED solver ----------
get_idx_tpl = lambda x: tuple(next(iter(x))[0][0][1])
op_indices = map(get_idx_tpl, util.get_fundamental_operators(h_imp))

index_converter = {}

for spin, orb in op_indices:
    # Bath degrees of freedom
    if isinstance(orb, str) and 'b_' in orb:
        index_converter[(spin, orb)] = ("bath", int(orb.split('_')[1]), "up" if spin == "up" else "down")
    # Local degrees of freedom
    else:
        index_converter[(spin, orb)] = ("loc", orb, "up" if spin == "up" else "down")

ed = PomerolED(index_converter, verbose = True)

# --------- Calculate the single-particle Green function ----------
ed.diagonalize(h_imp)
G_iw = ed.G_iw(dict([(v[0],v[1]) for v in gf_struct]), beta, n_iw)

# -------- Save in archive ---------
with HDFArchive("../results/pomerol.h5",'w') as res:
    res["G"] = G_iw
