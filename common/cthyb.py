execfile('../system.py')

from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from pytriqs.cthyb import Solver, version

# --------- Construct the CTHYB solver ----------
S = Solver(beta = beta, 
            gf_struct = dict(gf_struct),
            n_iw = n_iw,  
            n_tau = 100001)

# --------- Initialize G0_iw ----------

# FIXME Set tails explicitly
S.G0_iw['up'] << inverse(iOmega_n + mu + h)
S.G0_iw['dn'] << inverse(iOmega_n + mu - h)

for iw in Delta['up'].mesh:
    S.G0_iw['up'][iw] = 1.0 / (iw + mu + h - Delta['up'][iw])
    S.G0_iw['dn'][iw] = 1.0 / (iw + mu - h - Delta['dn'][iw])

# --------- Solve! ----------
S.solve(h_int=h_int,
        n_warmup_cycles = 1000,
        n_cycles = 10000000,
        length_cycle = 200,
       )

# -------- Save in archive ---------
if mpi.is_master_node():
    with HDFArchive("../results/cthyb.h5",'w') as res:
        res["G"] = S.G_iw
        res.create_group("Solver_Info")
        info_grp = res["Solver_Info"]
        info_grp["solver_name"] = "cthyb"
        info_grp["version"] = version.version
        info_grp["git_hash"] = version.cthyb_hash
        info_grp["triqs_git_hash"] = version.triqs_hash
