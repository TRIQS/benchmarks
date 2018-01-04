import sys, os
sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/../../common")
from model import *

from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

# --------- Construct the ED solver ----------
ed = TriqsExactDiagonalization(h_imp, fundamental_operators, beta)

# --------- Calculate the single-particle Green function ----------
G_iw = G0_iw.copy() 
for bl, idx_lst in gf_struct:
    for i, j in product(idx_lst, idx_lst):
        ed.set_g2_iwn(G_iw[bl][i,j], c(bl,i), c_dag(bl,j))

# -------- Save in archive ---------
with HDFArchive("../results/ed.h5",'w') as res:
    res["G"] = G_iw
