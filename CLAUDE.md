# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

TRIQS Solver Benchmarks: Systematic tests comparing quantum impurity solvers using the TRIQS library framework. Each model directory contains a complete benchmark for a specific impurity problem.

## Running Solver Scripts

Individual scripts:
```bash
cd Hubbard_Atom/scripts
mpirun -np 4 python cthyb                 # default: G + static observables
python atomdiag --measure G_w              # with additional observables
python ctint --measure chi3                # with chi3 measurement
```

All benchmarks via runner:
```bash
python run_benchmarks.py                          # all models
python run_benchmarks.py Hubbard_Atom             # one model
python run_benchmarks.py Hubbard_Atom ctint       # one model, one solver
python run_benchmarks.py --solver ctint           # one solver across all models
python run_benchmarks.py --dry-run                # print commands only
```

Common scripts in `common/` are symlinked into each model's `scripts/` directory.

## Architecture

### Model Directories
Each model contains:
- `model.py` - Defines physical parameters: `beta`, `h_int`, `G0_iw`, `gf_struct`, `n_iw`, and optionally hamiltonians (`h_0`, `h_imp`, `h_bath`, `h_coup`, `h_tot`), hybridization (`Delta_iw`), DLR mesh quantities (`G0_dlr_iw`, `Delta_dlr`, `dlr_iw_mesh`), real-frequency mesh (`w_mesh`, `broadening`), and dynamic interactions (`D0_func`, `Jperp_func`). Multi-orbital models also export `c_dag_vec`/`c_vec` operator vectors and `block_names`. The module-level docstring holds the markdown/LaTeX model description used as the report header.
- `scripts/` - Solver-specific scripts (symlinks to `common/` or custom)
- `results/` - HDF5 archives with solver output plus the generated `report.md`, `figures/`, and `tables/` artifacts (ignored on dev branches via .gitignore)

Models: Hubbard_Atom, SIAM_Discrete_Bath, SIAM_SemiCircular, SIAM_DynU, SIAM_Jperp, Dimer, Dimer_nn, Dimer_SOC, Trimer, Sr2RuO4, Sr2RuO4_SOC, La2CuO4, Plaquette, Plaquette_SemiCircular

### Solver Scripts Pattern
All solver scripts follow the same structure:
1. Import model parameters from `../model.py` via `from model import *`
2. Parse `--measure` arguments via `parse_measure_args(SUPPORTED)`
3. Auto-detect dynamic interactions (`D0_func`, `Jperp_func`) from model.py
4. Configure and run the solver
5. Collect observables (G, static_obs, optionally G_w/chi3/chi4)
6. Save via `save_results()` to `../results/SOLVER.h5`

### HDF5 Archive Structure
Results are stored in `results/SOLVER.h5` with standardized keys:

**Single-particle:** `G` (MeshImFreq/MeshDLRImFreq), `G_w` (MeshReFreq), `Sigma_w` (MeshReFreq)

**Two-particle (physical channel basis d/m/s/t):**
- `chi2_d`, `chi2_m`, `chi2_s`, `chi2_t` -- susceptibilities (bosonic Omega)
- `chi3_d`, `chi3_m`, `chi3_s`, `chi3_t` -- three-point correlators (Omega, nu)
- `chi4_d`, `chi4_m`, `chi4_s`, `chi4_t` -- connected G2c (Omega, nu, nu')

**Static and metadata:**
- `static_obs/density` -- n_a per block/orbital
- `static_obs/nn_ab` -- <n_a n_b> for all pairs
- `Solver_Info/` -- solver_name, solver_version, solver_git_hash, triqs_version, constr_params, solve_params, script, num_threads, run_time

### Common Utilities
- `common/save_utils.py` - `save_results()` standardized HDF5 output, `parse_measure_args()` for --measure CLI
- `common/channels.py` - Channel translations (pp/ph/xph) and spin decomposition (d/m/s/t)
- `common/analysis.py` - Shared analysis/plotting library: `load_all_results()`, `compute_sigma()`, `deviation_table()`, `deviation_vs_reference()`, `plot_iw_comparison()`, `plot_iw_scalar_comparison()`, `plot_w_comparison()`, `plot_static_obs_table()`, `plot_chi_contour()`. The `return_markdown=True` flag on table helpers yields a markdown string instead of printing.
- `common/build_report.py` - Per-model report generator. Produces `MODEL/results/report.md`, `figures/*.png`, `tables/*.md`.
- `common/build_summary.py` - Cross-model summary generator. Produces `summary/correctness.md`, `summary/coverage.md`, `summary/README.md`.
- `common/util.py` - `get_fundamental_operators(op)`: extracts annihilation operators from many-body operator expressions
- `common/util_mpi.py` - `mpi_print`: prints only on MPI master node
- `common/script_template` - Template for adding new solvers

### Common Solver Scripts
Not all models symlink all scripts; each model only includes the solvers applicable to it.
- `common/cthyb` - CT-HYB Monte Carlo (triqs_cthyb). Auto-detects `Jperp_func` for delta_interface.
- `common/ctseg` - CT-HYB segment picture (triqs_ctseg), Delta_tau interface with eigenbasis rotation. Auto-detects `D0_func`/`Jperp_func`.
- `common/ctint` - CT-INT Monte Carlo (triqs_ctint). Auto-detects `D0_func`/`Jperp_func`.
- `common/pyed` - Exact diagonalization. Supports `--measure chi3 chi4`.
- `common/pomerol` - Full exact diagonalization (pomerol2triqs)
- `common/atomdiag` - Atomic diagonalization (requires `h_tot`, `n_orb_bath`). Supports `--measure G_w`.
- `common/edipack` - EDIpack exact diagonalization (edipack2triqs, requires `h_tot`, `n_orb_bath`). Supports `--measure G_w chi2`. Provides G(iw), G(w), chi2_d (density), chi2_m (spin).
- `common/nrg` - NRG via nrgljubljana_interface (requires `w_mesh`)
- `common/alps_cthyb` - ALPS CT-HYB (via DCore)
- `common/w2dyn_cthyb` - w2dynamics CT-HYB
- `common/forktps` - Fork TPS solver (real-frequency)

### Custom Model Scripts
- `Hubbard_Atom/scripts/exact` - Analytical solution for Hubbard atom
- `Sr2RuO4/scripts/cthyb_truncated` - CT-HYB with full Hilbert space
- `Sr2RuO4/scripts/cthyb_truncation_benchmark` - Truncation level benchmarks
- `Plaquette_SemiCircular/scripts/cthyb_truncation_benchmark` - Truncation benchmarks

### Dynamic Interactions
Models with dynamic interactions export descriptors in model.py:
- `SIAM_DynU`: `D0_func` (density-density), `D0_block_pairs` -- auto-detected by ctint/ctseg
- `SIAM_Jperp`: `Jperp_func` (spin-spin) -- auto-detected by cthyb/ctint/ctseg

### Run Automation
- `benchmark_config.yaml` - Defines which solvers run on which models, with optional `--measure` args and a `reference_solver` per model for cross-model correctness comparisons
- `run_benchmarks.py` - Reads config and runs all benchmarks, captures logs to `MODEL/results/SOLVER.log`, produces `benchmark_results.json`. Runs each solver from `MODEL/scripts/` with a 1-hour timeout.

Config format:
```yaml
defaults:
  n_mpi_ranks: 4
  reference_solver_priority: [exact, pyed, pomerol, edipack, atomdiag]  # fallback if per-model reference_solver is absent
models:
  ModelName:
    reference_solver: pyed                        # used by build_summary.py for deviations
    runs:
      - solver_name                              # string form (no measure args)
      - {solver: solver_name, measure: all}      # dict with measure (string or list)
      - {solver: solver_name, measure: [chi2, chi3]}
    n_mpi_ranks: 16  # optional per-model override
```

### Report Generation
- `python common/build_report.py [MODEL ...]` - Build `MODEL/results/report.md` with embedded plots and markdown deviation tables. Imports each model's `model.py` with cwd set to `MODEL/` (to satisfy model.py's `sys.path.append(os.getcwd() + '/../common')`). Requires `G0_iw` on the module for Sigma derivation.
- `python common/build_summary.py` - Build the cross-model `summary/` directory. Uses `reference_solver` from the config for max-norm deviation columns. Models without a reference appear with an `_all_` placeholder row.

## Adding a New Solver

1. Copy `common/script_template` to `common/your_solver`
2. Implement solver initialization, solve call, and observable collection
3. Use `save_results()` for standardized output
4. Symlink into applicable model directories: `ln -s ../../common/your_solver scripts/your_solver`
5. Add to `benchmark_config.yaml`

## Adding a New Model

1. Create a new directory with `model.py` defining `beta`, `h_int`, `G0_iw`, `gf_struct`, `n_iw`. Add a module docstring at the top with a markdown/LaTeX description — it will be the header of the generated report.
2. Create `scripts/` and `results/` subdirectories
3. Symlink common solver scripts: `ln -s ../../common/ctint scripts/ctint`
4. Add to `benchmark_config.yaml` (with a `reference_solver` if you have an ED option available)
5. After running the benchmarks, regenerate artifacts: `python common/build_report.py NewModel` and `python common/build_summary.py`

## Conventions and Gotchas

- **Symlinks only**: Always symlink common solvers into `MODEL/scripts/`; never hardcopy. Centralized bugfixes must propagate.
- **Path setup in scripts**: All solver scripts do `sys.path.append(os.getcwd() + '/..')` to import `model.py` from the parent directory. Scripts must be run from `MODEL/scripts/`.
- **save_results() is MPI-aware**: Only saves on master node. Non-master ranks silently return. It also captures the full script source via `inspect.getsource(__main__)`.
- **DLR requirement for ctint**: ctint requires `G0_dlr_iw` (not `G0_iw`) for solver initialization.
- **Eigenbasis rotation for ctseg**: Unlike cthyb, ctseg requires rotating G0 into the eigenbasis of h0 (extracted via hermitian tail fitting) before Fourier transform to tau.
- **Static obs key naming**: Single-orbital models use short keys ("up", "dn"); multi-orbital use indexed keys ("up_0", "up_1"). Conditional: `key = f"{bl}_{i}" if bl_size > 1 else bl`.
- **`--measure all`**: Expands to full SUPPORTED list within `parse_measure_args()`; solvers never see the literal "all" keyword.
- **Chi3 DLR2D in ctint**: ctint stores chi3 in compressed DLR2D NFFT form; must convert via `make_gf_imfreq(make_gf_dlr2d(...))` for analysis.
- **Real-frequency solvers**: atomdiag/edipack/nrg/forktps only compute G_w if `model.py` defines `w_mesh`, `broadening`, and `w_window`.
- **Channel inter-translations**: `channels.py` pp/ph/xph translations raise `NotImplementedError`; only spin decomposition (d/m/s/t) is implemented.
- **DLR2D branch**: triqs/cthyb/ctseg/ctint now build from `DLR2D` branch. cthyb/ctseg auto-determine `n_warmup_cycles` and `length_cycle` -- do not set them manually.
- **Docker usage**: Build with `docker build -t solver_benchmarks .` (use `--build-arg NCORES=16` for faster builds). For short tasks (single solver, `--dry-run`, `build_report.py`, `build_summary.py`) run foreground: `docker run --rm -u 0:0 -v $(pwd):/home/triqs/benchmarks solver_benchmarks <cmd>`. **For a full `run_benchmarks.py` pass (hours long) always run detached** so the container survives terminal/SSH disconnect: `docker run -d --name solver_bench_run -u 0:0 -v $(pwd):/home/triqs/benchmarks solver_benchmarks bash -c "python -u run_benchmarks.py 2>&1 && python -u common/build_report.py && python -u common/build_summary.py"`. A foreground `docker run --rm` for the full pass was observed to be killed when its client/parent-shell went away, leaving the container orphaned and auto-removed. `python -u` is required for unbuffered `docker logs` progress. Follow with `docker logs -f solver_bench_run`, clean up afterwards with `docker rm solver_bench_run`. The `-u 0:0` is needed due to Docker user namespace remapping causing UID mismatch with bind mounts. Results are written to `MODEL/results/` on the host via the bind mount. For interactive use: `docker run --rm -it -u 0:0 -v $(pwd):/home/triqs/benchmarks solver_benchmarks`.
- **DCore MPI conflict**: `save_utils.py` defers `from triqs.utility import mpi` into function body to avoid polluting `sys.modules` before DCore's `raise_if_mpi_imported()` check.
