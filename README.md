TRIQS Bench
===========

TRIQS Bench provides a variety of benchmark systems for quantum impurity solvers.
We use the [TRIQS library](https://www.triqs.ipht.fr/) as a framework for defining the quantum systems
and for the comparison of different results.

Each directory fully defines one specific impurity model and contains the following:

* **README.rst** - A description of the impurity model
* **system.py** - Define the local Hamiltonian and the hybridization function of the impurity model
* **scripts** - Contains one script for each impurity solver (e.g. CTHYB.py)
