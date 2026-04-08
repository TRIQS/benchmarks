TRIQS Solver Benchmarks
=======================

This repository provides systematic tests and benchmarks of various quantum impurity solvers
using the Python interface of the [TRIQS library](https://triqs.github.io/triqs) as a framework.

Each directory defines one specific impurity model and contains:

* **model.py** -- Hamiltonian, hybridization function, and Green function structure.
* **scripts/** -- One script per applicable impurity solver (symlinks to `common/`).
* **results/** -- HDF5 archives with solver output.
* **notebook.ipynb** -- Jupyter notebook with model description and comparison of results.

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

Models
------

| Model | Description |
|---|---|
| **Hubbard_Atom** | Single site with Coulomb repulsion, chemical potential, and Zeeman field |
| **SIAM_Discrete_Bath** | Single impurity Anderson model coupled to discrete bath levels |
| **SIAM_Wide_Band** | Single impurity Anderson model with semicircular (wide-band) hybridization |
| **SIAM_DynU** | SIAM with dynamic density-density interaction $D_0(i\omega)$ |
| **SIAM_Jperp** | SIAM with dynamic spin-spin interaction $J_\perp(i\omega)$ |
| **Dimer** | Two-orbital Kanamori impurity coupled to discrete bath |
| **Dimer_nn** | Two-orbital density-density impurity coupled to discrete bath |
| **Dimer_SOC** | Two-orbital impurity with spin-orbit coupling (single spin-orbital block) |
| **Trimer** | Three-orbital Kanamori impurity coupled to discrete bath |
| **Plaquette** | Four-site cluster with Kanamori interaction and discrete bath |
| **Plaquette_Wide_Band** | Four-site cluster with wide-band hybridization |
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

Adding Your Solver
------------------

1. Copy `common/script_template` to `common/your_solver`.
2. Implement solver initialization, solve call, and observable collection.
3. Use `save_results()` from `common/save_utils.py` for standardized HDF5 output.
4. Symlink into applicable model directories: `ln -s ../../common/your_solver scripts/your_solver`
5. Add entries to `benchmark_config.yaml`.

For questions, [open an issue](https://github.com/TRIQS/benchmarks/issues) or contact [Nils Wentzell](mailto:nils.wentzell@gmail.com).
