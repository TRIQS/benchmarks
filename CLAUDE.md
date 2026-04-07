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
- `model.py` - Defines physical parameters: `beta`, `h_int`, `G0_iw`, `gf_struct`, `n_iw`, and optionally hamiltonians (`h_0`, `h_imp`, `h_bath`, `h_coup`, `h_tot`), hybridization (`Delta_iw`), DLR mesh quantities (`G0_dlr_iw`, `Delta_dlr`, `dlr_iw_mesh`), real-frequency mesh (`w_mesh`, `broadening`), and dynamic interactions (`D0_func`, `Jperp_func`). Multi-orbital models also export `c_dag_vec`/`c_vec` operator vectors and `block_names`.
- `scripts/` - Solver-specific scripts (symlinks to `common/` or custom)
- `results/` - HDF5 archives with solver output (ignored on dev branches via .gitignore)
- `notebook.ipynb` + `notebook.py` - Analysis and comparison (jupytext paired)

Models: Hubbard_Atom, SIAM_Discrete_Bath, SIAM_Wide_Band, SIAM_DynU, SIAM_Jperp, Dimer, Dimer_nn, Dimer_SOC, Trimer, Sr2RuO4, Sr2RuO4_SOC, La2CuO4, Plaquette, Plaquette_Wide_Band

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
- `common/analysis.py` - Shared analysis/plotting library: `load_all_results()`, `compute_sigma()`, `deviation_table()`, `plot_iw_comparison()`, `plot_w_comparison()`, `plot_static_obs_table()`, `plot_chi_contour()`
- `common/util.py` - `get_fundamental_operators(op)`: extracts annihilation operators from many-body operator expressions
- `common/util_mpi.py` - `mpi_print`: prints only on MPI master node
- `common/plot.py` - Legacy plotting helpers (kept for backward compatibility)
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
- `Plaquette_Wide_Band/scripts/cthyb_truncation_benchmark` - Truncation benchmarks

### Dynamic Interactions
Models with dynamic interactions export descriptors in model.py:
- `SIAM_DynU`: `D0_func` (density-density), `D0_block_pairs` -- auto-detected by ctint/ctseg
- `SIAM_Jperp`: `Jperp_func` (spin-spin) -- auto-detected by cthyb/ctint/ctseg

### Run Automation
- `benchmark_config.yaml` - Defines which solvers run on which models, with optional `--measure` args
- `run_benchmarks.py` - Reads config and runs all benchmarks, captures logs, produces summary JSON

## Adding a New Solver

1. Copy `common/script_template` to `common/your_solver`
2. Implement solver initialization, solve call, and observable collection
3. Use `save_results()` for standardized output
4. Symlink into applicable model directories: `ln -s ../../common/your_solver scripts/your_solver`
5. Add to `benchmark_config.yaml`

## Adding a New Model

1. Create a new directory with `model.py` defining `beta`, `h_int`, `G0_iw`, `gf_struct`, `n_iw`
2. Create `scripts/` and `results/` subdirectories
3. Symlink common solver scripts: `ln -s ../../common/ctint scripts/ctint`
4. Create `notebook.py` from template, convert with `jupytext --to ipynb notebook.py`
5. Add to `benchmark_config.yaml`

## Impurity Solvers Used
- triqs_cthyb - CT-HYB Monte Carlo
- triqs_ctseg - CT-HYB segment picture
- triqs_ctint - CT-INT Monte Carlo
- pyed - Exact diagonalization
- pomerol - Full ED
- edipack2triqs - EDIpack ED
- nrgljubljana - NRG
- w2dynamics - w2dyn CT-HYB
- alps_cthyb - ALPS CT-HYB
- forktps - Fork TPS
