TRIQS Solver benchmarks
=======================

This repository provides systematic tests and benchmarks of various impurity solvers
for a set of reference impurity models. We use the Python interface of the
[TRIQS library](https://triqs.github.io/triqs) as a framework.

Each directory fully defines one specific impurity model and contains the following:

* **notebook.ipynb** - IPython Notebook with a description of the impurity model and analysis of results.
* **model.py** - Defines the local Hamiltonian and the hybridization function of the impurity model.
* **scripts** - Contains one script for each applicable impurity solver (e.g. cthyb, pyed, ...).
* **results** - Contains one hdf5 archive for each impurity solver (e.g. cthyb.h5).
		They consist of the results for the Green function and additional solver information.

Models
------

* [**Hubbard_Atom**](https://github.com/TRIQS/benchmarks/blob/master/Hubbard_Atom/notebook.ipynb) A single atomic level with a Coulomb repulsion, a chemical potential and a Zeeman splitting term
* [**SIAM_Discrete_Bath**](https://github.com/TRIQS/benchmarks/blob/master/SIAM_Discrete_Bath/notebook.ipynb) A dimer with spin-orbit coupling and density-density interaction coupled to two discrete bath states
* [**SIAM_Wide_Band**](https://github.com/TRIQS/benchmarks/blob/master/SIAM_Wide_Band/notebook.ipynb) A dimer with spin-orbit coupling and density-density interaction coupled to two discrete bath states
* [**Dimer**](https://github.com/TRIQS/benchmarks/blob/master/Dimer/notebook.ipynb) A dimer with Kanamori-Interaction coupled to two discrete bath states
* [**Dimer_SOC**](https://github.com/TRIQS/benchmarks/blob/master/Dimer_SOC/notebook.ipynb) A dimer with spin-orbit coupling and density-density interaction coupled to two discrete bath states
* [**Trimer**](https://github.com/TRIQS/benchmarks/blob/master/Trimer/notebook.ipynb) A trimer with Kanamori-Interaction coupled to three discrete bath states
* [**Sr2RuO4**](https://github.com/TRIQS/benchmarks/blob/master/Sr2RuO4/notebook.ipynb) An effective 3-band impurity model for Sr2RuO4
* [**Sr2RuO4_SOC**](https://github.com/TRIQS/benchmarks/blob/master/Sr2RuO4/notebook.ipynb) An effective 3-band impurity model for Sr2RuO4 including spin-orbit coupling

Impurity Solvers
----------------

* [**triqs_cthyb**](https://triqs.github.io/cthyb) - Continuous-time hybridization-expansion quantum Monte-Carlo code based on TRIQS.<br/>  Maintainer: [Nils Wentzell](mailto:nwentzell@flatironinstitute.com)
* [**triqs_ctseg**](https://triqs.github.io/ctseg) (private) - Continuous-time hybridization-expansion quantum Monte-Carlo code in the segment picture.<br/>  Maintainer: [Thomas Ayral](mailto:th.ayral@gmail.com)
* [**triqs_ctint**](https://triqs.github.io/ctint) (private) - Continuous-time interaction-expansion quantum Monte-Carlo code based on TRIQS.<br/>  Maintainer: [Nils Wentzell](mailto:nwentzell@flatironinstitute.com)
* [**pyed**](https://github.com/hugostrand/pyed) - Exact diagonalization solver for finite quantum systems based on TRIQS. <br/>  Maintainer: [Hugo Strand](mailto:hstrand@flatironinstitute.org)
* [**pomerol**](https://github.com/aeantipov/pomerol) - An exact diagonalization (full-ED) code written in C++ aimed at solving condensed matter second-quantized models of interacting fermions on finite size lattices at finite temperatures. It is designed to produce single and two-particle Greens functions. ([TRIQS Interface](https://github.com/krivenko/pomerol2triqs)).<br/>  Maintainer: [Andrey Antipov](mailto:andrey.e.antipov@gmail.com)
* [**w2dynamics**](https://github.com/TRIQS/w2dynamics/w2dynamics) - A continuous-time hybridization expansion impurity solver contained in the w2dynamics software package ([TRIQS interface](https://triqs.github.io/w2dynamics_interface)).<br/>  Maintainer: [Andreas Hausoel](mailto:andreas.hausoel@physik.uni-wuerzburg.de)

Adding your Solver
------------------

If you would like to add your impurity solver to this benchmark project, please start from the instructions given in the [script_template](https://github.com/TRIQS/benchmarks/blob/master/common/script_template) file.
For questions feel free to [contact me](mailto:nils.wentzell@gmail.com) or post an issue.
