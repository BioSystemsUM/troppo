|License| |PyPI version|

TROPPO
============

*Troppo* (Tissue-specific RecOnstruction and Phenotype Prediction using Omics data) is a Python package containing methods
for tissue specific reconstruction to use with constraint-based models. The main purpose of this package is to provide
an open-source framework which is both modular and flexible to be integrated with other packages, such as cobrapy, framed
or cameo whose already provide generic data structures to interpret metabolic models.

A (MI)LP solver is required to use most of the present methods. The current methods support optlang, which in turn allow
the use of solvers like CPLEX or Gurobi.

The current methods implemented are:
    - FastCORE
    - CORDA
    - GIMME
    - (t)INIT
    - iMAT

Methods to be implemented later:
    - MBA
    - mCADRE
    - PRIME

Documentation
~~~~~~~~~~~~~



Instalation from PyPI (stable releases)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    pip install troppo

Instalation from github (latest development release)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    pip install https://github.com/BioSystemsUM/troppo

Credits and License
~~~~~~~~~~~~~~~~~~~

Developed at the Centre of Biological Engineering, University of Minho

Released under the GNU Public License (version 3.0).


.. |License| image:: https://img.shields.io/badge/license-GPL%20v3.0-blue.svg
   :target: https://opensource.org/licenses/GPL-3.0
.. |PyPI version| image:: https://badge.fury.io/py/troppo.svg
   :target: https://badge.fury.io/py/troppo