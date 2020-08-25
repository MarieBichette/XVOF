# XVOF
One dimensional hydro code for testing xfem method. 
This software is an experimental work trying to model spall fracture using enrichment and partition of unity 
method in a finite volume code. Hansbo & Hansbo enrichment is used.

## Code description

#### Available physics

###### Bulk behavior
Material behavior includes simple models for hydrodynamics, elasticity, plasticity.
- For hydrodynamics, the equation of state is a MieGrüneisen type eos
- Linear elasticity and perfect plasticity is implemented
- **No** viscous plasticity is considered

###### Failure models
Two rupture models are available :

1) *imposed pressure* : once the cell reach rupture criterion, the pressure is imposed at the value specified in data file

2) *enrichment* : once the cell reach rupture criterion, kinematic and thermodynamic fields are enriched according to the Hansbo \& Hansbo idea.
Mass matrix can be lumped (Menouillard or total lumping) to improve calculation time.
Penalty method is used to avoid the interpenetration of the discontinuity boundaries.

Damage can be modeled with a cohesive zone model with specific cohesive law. But no bulk damage is implemented

#### Non regression procedure
This code is versioned using GitHub. A non regression procedure is launched  by Travis at each *push* on the GitHub repository 
* A unittest collection
* Integration tests for hydrodynamics, elasticity, plasticity with and without enrichment (result should be compared to reference.hdf5)

Coverage is also tested within the non regression procedure.

#### Future work
* Check CFL stability with respect to mass matrix.
* Implement other contact treatments.
* Implement simple damage models

#### Miscellaneous
This code is now ported in Python3.7, using type hints for local variables and method arguments.

## Installation
- Download the GitHub XVOF repository
- To install the lib C for equation of state computation :
from XVOF repository type : ***(python3.7 -m) pip install -e .***

## Test case creation
Each case is composed of :
1) Input data, stored in the file XDATA.json.
2) A mesh file, which is a .txt file containing the coordinates of nodes.

To launch a test case :
- Create a data set "XDATA.json" and a meshfile "mesh.txt"
- Type : ***OMP_NUM_THREADS=2 python3.7 XtendedFiniteVolume <case-repository>***
