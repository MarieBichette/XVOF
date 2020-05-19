Case creation
============================
Each case is composed of :

1) Input data, stored in the file XDATA.json.
2) A mesh file, which is a .txt file containing the coordinates of nodes.

To launch a test case :

- Create a data set "XDATA.json" and a meshfile "mesh.txt"
- From the test repository, type :

.. code-block::

    OMP_NUM_THREADS=2 python XtendedFiniteVolume <case-repository>

By default, the external lib C is used if it has been previously installed.

To enforce the internal computation of the equation of state (with python module), type

.. code-block::

    python XtendedFiniteVolume <case-repository> --use-internal-solver