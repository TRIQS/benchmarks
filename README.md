TRIQS Bench
===========

TRIQS Bench provides a variety of benchmark models for quantum impurity solvers.
We use the [TRIQS library](https://www.triqs.ipht.fr/) as a framework for defining the quantum systems
and for the comparison of different results.

Each directory fully defines one specific impurity model and contains the following:

* **notebook.ipynb** - IPython Notebook with a description of the impurity model and analysis of results.
* **model.py** - Defines the local Hamiltonian and the hybridization function of the impurity model.
* **scripts** - Contains one script for each impurity solver (e.g. cthyb.py).
* **results** - Contains one hdf5 archive for each impurity solver (e.g. cthyb.h5).
		They consist of the results for the Green function and additional solver information.
