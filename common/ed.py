execfile('../system.py')

from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

# --------- Construct the ED solver ----------
ed = TriqsExactDiagonalization(h_imp, fundamental_operators, beta)

# --------- Calculate the single-particle Green function ----------
G_iw = Gf_from_struct(mesh=iw_mesh, struct=gf_struct)
G_iw['up'] << inverse(iOmega_n + mu + h)
G_iw['dn'] << inverse(iOmega_n + mu - h)

for bl, idx_lst in gf_struct:
    for i, j in product(idx_lst, idx_lst):
        ed.set_g2_iwn(G_iw[bl][i,j], c(bl,i), c_dag(bl,j))

# -------- Save in archive ---------
with HDFArchive("../results/ed.h5",'w') as res:
    res["G"] = G_iw
