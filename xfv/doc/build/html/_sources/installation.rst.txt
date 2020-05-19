Installation
============================
- Download the GitHub XVOF repository (`git clone` or download)
- Install pyenv with a python version > 3.7 (follow instructions in https://github.com/pyenv/pyenv#installation). This will enable to modify the python of pyenv instead of the system python
- Install cmake and swig :

.. code-block::

   apt-get install cmake
   apt-get install swig

Cmake will enable to build the C library for the eos computation and SWIG translates it as a python module
Cmake version should be > 3.12.

- To activate python of pyenv :

from XVOF repository type

.. code-block::

   pyenv local 3.7.7

(python3.8 can currently create bug when installing the library C for equation of state, so it is recommended to use 3.7 version)

- To install the XVOF code and dependencies :

from XVOF repository type :

.. code-block::

   pip install --upgrade pip
   pip install -e .

This will also build the lib C for equation of state computation.