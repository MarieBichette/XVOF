# XVOF
One dimensional hydro code for testing xfem method. 
This software is an experimental work trying to model spall fracture using enrichment and partition of unity 
method in a finite volume code. Hansbo & Hansbo enrichment is used.

## Available physics
Material behavior includes simple models for hydrodynamics, elasticity, plasticity.

Two rupture models are available :

1) *imposed pressure* : once the cell reach rupture criterion, the pressure is imposed at the value specified in data file

2) *enrichment* : once the cell reach rupture criterion, kinematic and thermodynamic fields are enriched according to the Hansbo \& Hansbo idea.
Mass matrix can be lumped (Menouillard or total lumping) to improve calculation time.

Damage can be modeled with a cohesive zone model with specific cohesive law.

Penalty method is used to avoid the interpenetration of the discontinuity boundaries.

## Case creation
Each case is composed of :

1) Input data, stored in the file XDATA.json.
2) A mesh file, which is a .txt file containing the coordinates of nodes.

## Non regression procedure
* A unittest collection
* Integration tests for hydrodynamics, elasticity, plasticity with and without enrichment (result should be compared to reference.hdf5)

The complete procedure is executed by Travis at each *push* of the Git repository. Coverage is also tested within the non regression procedure.

## Future work
* Check CFL stability with respect to mass matrix.
* Use a C module to call the equation of state function in order to spare computation time.
* Implement other contact treatments.

## Miscellaneous
This code is now ported in Python3.7, using type hints for local variables and method arguments.

To launch a profiling under ipython type:
%run -p -r -s module -s cumulative XtentedFiniteVolume.py
