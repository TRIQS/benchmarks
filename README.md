TRIQS Solver Benchmarks
=======================

This repository provides systematic tests and benchmarks of various quantum impurity solvers
using the Python interface of the [TRIQS library](https://triqs.github.io/triqs) as a framework.

Each directory defines one specific impurity model and contains:

* **model.py** -- Hamiltonian, hybridization function, and Green function structure. The module docstring carries a markdown/LaTeX description of the model used as the header of the generated report.
* **scripts/** -- One script per applicable impurity solver (symlinks to `common/`).
* **results/** -- HDF5 archives with solver output plus a generated `report.md` containing per-observable plots and deviation tables.

Running Benchmarks
------------------

Individual solver scripts:
```bash
cd Hubbard_Atom/scripts
mpirun -np 4 python cthyb                 # G(iw) + static observables
python atomdiag --measure G_w              # with real-frequency measurement
python ctint --measure chi2 chi3           # with two-particle measurements
```

All benchmarks via the runner:
```bash
python run_benchmarks.py                          # all models
python run_benchmarks.py Hubbard_Atom             # one model
python run_benchmarks.py Hubbard_Atom ctint       # one model, one solver
python run_benchmarks.py --solver ctint           # one solver across all models
python run_benchmarks.py --dry-run                # print commands only
python run_benchmarks.py --slurm                  # generate SLURM job scripts
```

The runner reads `benchmark_config.yaml` for the list of models, solvers, and measurement flags.

Analysis & Reports
------------------

Comparison plots and deviation tables are produced as static markdown + PNG
artifacts under each model's `results/` directory. No Jupyter kernel required
-- GitHub renders the LaTeX and images inline.

Per-model reports (one markdown file per model, with plots and deviation tables
for every observable):
```bash
python common/build_report.py                     # all models with results
python common/build_report.py Hubbard_Atom Dimer  # selected models
```

This writes `MODEL/results/report.md`, `MODEL/results/figures/*.png`, and
`MODEL/results/tables/*.md` for each model.

Cross-model summary (correctness and coverage matrices across all models and
solvers):
```bash
python common/build_summary.py
```

This writes `summary/correctness.md`, `summary/coverage.md`, and
`summary/README.md`. The reference solver per model is read from the
`reference_solver` key in `benchmark_config.yaml`; if absent, the runner falls
back to the `reference_solver_priority` list in the `defaults:` section.

Models
------

| Model | Description |
|---|---|
| **Hubbard_Atom** | Single site with Coulomb repulsion, chemical potential, and Zeeman field |
| **SIAM_Discrete_Bath** | Single impurity Anderson model coupled to discrete bath levels |
| **SIAM_SemiCircular** | Single impurity Anderson model with semicircular hybridization |
| **SIAM_DynU** | SIAM with dynamic density-density interaction $D_0(i\omega)$ |
| **SIAM_Jperp** | SIAM with dynamic spin-spin interaction $J_\perp(i\omega)$ |
| **Dimer** | Two-orbital Kanamori impurity coupled to discrete bath |
| **Dimer_nn** | Two-orbital density-density impurity coupled to discrete bath |
| **Dimer_SOC** | Two-orbital impurity with spin-orbit coupling (single spin-orbital block) |
| **Trimer** | Three-orbital Kanamori impurity coupled to discrete bath |
| **Plaquette** | Four-site cluster with Kanamori interaction and discrete bath |
| **Plaquette_SemiCircular** | Four-site cluster with semicircular hybridization |
| **Sr2RuO4** | Three-band model for Sr$_2$RuO$_4$ from Wannier90 |
| **Sr2RuO4_SOC** | Three-band Sr$_2$RuO$_4$ with spin-orbit coupling |
| **La2CuO4** | Single-band model for La$_2$CuO$_4$ from Wannier90 |

Impurity Solvers
----------------

| Solver | Script | Type |
|---|---|---|
| [triqs_cthyb](https://triqs.github.io/cthyb) | `cthyb` | CT-HYB Monte Carlo |
| [triqs_ctseg](https://triqs.github.io/ctseg) | `ctseg` | CT-HYB segment picture |
| [triqs_ctint](https://triqs.github.io/ctint) | `ctint` | CT-INT Monte Carlo |
| [pyed](https://github.com/hugostrand/pyed) | `pyed` | Exact diagonalization |
| [pomerol](https://github.com/krivenko/pomerol2triqs) | `pomerol` | Full ED (two-particle) |
| [edipack2triqs](https://github.com/EDIpack/edipack2triqs) | `edipack` | EDIpack ED |
| [nrgljubljana_interface](https://github.com/TRIQS/nrgljubljana_interface) | `nrg` | NRG |
| [w2dynamics](https://triqs.github.io/w2dynamics_interface) | `w2dyn_cthyb` | w2dyn CT-HYB |
| [ALPS/CT-HYB](https://github.com/ALPSCore/CT-HYB) | `alps_cthyb` | ALPS CT-HYB (via DCore) |
| [ForkTPS](https://github.com/TRIQS/forktps) | `forktps` | Fork TPS (real-frequency) |
| triqs.atom_diag | `atomdiag` | Atomic ED |

Docker
------

A Dockerfile is provided that builds all supported solvers from source. This
gives a reproducible environment without needing to install TRIQS and the
solvers locally.

Build the image (takes ~30 min with 10 cores):
```bash
docker build -t solver_benchmarks .
docker build -t solver_benchmarks --build-arg NCORES=16 .   # more cores
```

Run benchmarks by bind-mounting the repository into the container:
```bash
docker run --rm -u 0:0 -v $(pwd):/home/triqs/benchmarks solver_benchmarks \
    python run_benchmarks.py                                # all benchmarks
docker run --rm -u 0:0 -v $(pwd):/home/triqs/benchmarks solver_benchmarks \
    python run_benchmarks.py Hubbard_Atom exact             # single solver
docker run --rm -u 0:0 -v $(pwd):/home/triqs/benchmarks solver_benchmarks \
    python run_benchmarks.py --dry-run                      # preview only
```

Results are written to `MODEL/results/` on the host via the bind mount.

> **Note:** `-u 0:0` is required because Docker user namespace remapping can
> cause a UID mismatch between the container user and the bind-mounted files.

For an interactive shell inside the container:
```bash
docker run --rm -it -u 0:0 -v $(pwd):/home/triqs/benchmarks solver_benchmarks
```

Adding Your Solver
------------------

1. Copy `common/script_template` to `common/your_solver`.
2. Implement solver initialization, solve call, and observable collection.
3. Use `save_results()` from `common/save_utils.py` for standardized HDF5 output.
4. Symlink into applicable model directories: `ln -s ../../common/your_solver scripts/your_solver`
5. Add entries to `benchmark_config.yaml`.
6. Regenerate reports so your solver appears in the comparison: `python common/build_report.py` (and `python common/build_summary.py` for the cross-model overview).

For questions, [open an issue](https://github.com/TRIQS/benchmarks/issues) or contact [Nils Wentzell](mailto:nils.wentzell@gmail.com).
