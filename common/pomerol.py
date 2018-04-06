import sys, os
sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/../../common")
from model import *
import util

from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from pomerol2triqs import PomerolED

# --------- Construct the ED solver ----------
index_converter = {}
# Local degrees of freedom
index_converter.update({(spin, orb) : ("loc", orb, "down" if spin == "dn" else "up") 
    for spin, orb in product(spin_names,orb_names)})
# Bath degrees of freedom
index_converter.update({(spin, 'b_' + str(orb)) : ("bath", orb, "down" if spin == "dn" else "up")
    for spin, orb in product(spin_names,orb_names)})

ed = PomerolED(index_converter, verbose = True)

# --------- Calculate the single-particle Green function ----------
ed.diagonalize(h_imp)
G_iw = ed.G_iw(dict([(v[0],v[1]) for v in gf_struct]), beta, n_iw)

# -------- Save in archive ---------
with HDFArchive("../results/pomerol.h5",'w') as res:
    res["G"] = G_iw
