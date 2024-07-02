Troppo
==============================

Troppo (Tissue-specific RecOnstruction and Phenotype Prediction using Omics data) is a Python package containing methods for tissue specific reconstruction to use with constraint-based models. The main purpose of this package is to provide an open-source framework which is both modular and flexible to be integrated with other packages, such as cobrapy, framed or cameo whose already provide generic data structures to interpret metabolic models.

A (MI)LP solver is required to use most of the present methods. The current methods support optlang, which in turn allow the use of solvers like CPLEX or Gurobi.

The current methods implemented are:

* FastCORE
* CORDA
* GIMME
* (t)INIT
* iMAT

.. toctree::
   :maxdepth: 3

   installation
   troppo.methods
   troppo.omics
   troppo.tasks
   troppo.utilities
   troppo.validation
