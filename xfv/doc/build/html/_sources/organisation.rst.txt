Code organisation
==========================
The XVOF repository contains a repository ``xfv/`` and different programs to be executed and parameter files to be used when modifications are pushed on the GitHub repository.

The ``xfv/`` is organized as follows :

- ``data/`` : available material data for constitutive model and an example of data set for case creation
- ``doc/`` : scripts to build code documentation (sphinx) and the html doc
- ``post_processing/`` : some scripts to plot the results (field evolution, space-time diagram, free surface velocity, ...)
- ``src/`` : source of the program
- ``tests/`` : non regression tests (unit tests data set and integration tests) and other tests
- ``generate_mesh.py`` : program to generate the mesh for a test case
- ``XtendedFiniteVolume.py`` : the executable to launch