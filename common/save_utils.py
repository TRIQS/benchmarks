import argparse
import inspect

from h5 import HDFArchive
from triqs.utility import mpi


def parse_measure_args(supported):
    """Parse --measure arguments from command line.

    Parameters
    ----------
    supported : list of str
        Observable names this solver can measure beyond defaults
        (e.g. ['chi3', 'chi4', 'G_w']).

    Returns
    -------
    argparse.Namespace with attribute `measure` (list of str).
    """
    parser = argparse.ArgumentParser(description="Run solver benchmark")
    parser.add_argument('--measure', nargs='+', default=[],
                        help="Additional observables to measure: " +
                             ", ".join(supported) + ", or 'all'")
    args = parser.parse_args()

    # Expand 'all' to the full supported list
    if 'all' in args.measure:
        args.measure = list(supported)
    else:
        for m in args.measure:
            if m not in supported:
                parser.error(f"Unknown observable '{m}'. "
                             f"Supported: {supported + ['all']}")

    return args


def save_results(filepath, solver_name, solver_version, solver_git_hash,
                 constr_params=None, solve_params=None, run_time=None,
                 **observables):
    """Save solver results to an HDF5 archive with standard metadata.

    Each entry in `observables` is stored as a top-level key.
    Nested dicts (e.g. static_obs) are stored as HDF5 subgroups.

    Parameters
    ----------
    filepath : str
        Output path, e.g. "../results/ctint.h5"
    solver_name : str
        e.g. "triqs_ctint"
    solver_version : str
    solver_git_hash : str
    constr_params : dict, optional
    solve_params : dict, optional
    run_time : float, optional
    **observables : dict
        Keys like G, G_w, chi3_d, chi3_m, static_obs, etc.
    """
    if not mpi.is_master_node():
        return

    import __main__

    with HDFArchive(filepath, 'w') as ar:
        # Write observables
        for key, val in observables.items():
            if isinstance(val, dict):
                ar.create_group(key)
                for subkey, subval in val.items():
                    ar[key][subkey] = subval
            else:
                ar[key] = val

        # Write metadata
        ar.create_group("Solver_Info")
        info = ar["Solver_Info"]
        info["solver_name"] = solver_name
        info["solver_version"] = solver_version
        info["solver_git_hash"] = solver_git_hash

        from triqs import version as triqs_version
        info["triqs_version"] = triqs_version.version
        info["triqs_git_hash"] = triqs_version.git_hash

        if constr_params is not None:
            info["constr_params"] = constr_params
        if solve_params is not None:
            info["solve_params"] = solve_params
        if run_time is not None:
            info["run_time"] = run_time

        info["script"] = inspect.getsource(__main__)
        info["num_threads"] = mpi.world.Get_size()
