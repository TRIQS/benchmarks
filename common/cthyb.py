execfile('../model.py')

from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from pytriqs.cthyb import Solver, version

# --------- Construct the CTHYB solver ----------
S = Solver(beta = beta, 
            gf_struct = dict(gf_struct),
            n_iw = n_iw,  
            n_tau = 100001)

# --------- Initialize G0_iw ----------
S.G0_iw << G0_iw

# --------- Solve! ----------
S.solve(h_int=h_int,
        n_warmup_cycles = 1000,
        n_cycles = 10000000,
        length_cycle = 200,
       )

# -------- Save in archive ---------
if mpi.is_master_node():
    with HDFArchive("../results/cthyb.h5",'w') as results:
        results["G"] = S.G_iw

        import inspect
        import __main__
        results.create_group("Solver_Info")
        info_grp = results["Solver_Info"]
        info_grp["solver_name"] = "triqs_ctint"
        info_grp["version"] = version.version
        info_grp["git_hash"] = version.cthyb_hash
        info_grp["triqs_git_hash"] = version.triqs_hash
        info_grp["script"] = inspect.getsource(__main__)
