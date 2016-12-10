# XVOF
One dimensional hydro code for testing xfem method

This software is an experimental work trying to modelize spall rupture thanks to enrichment and partition of unity 
method in a finite volume code.

Input data are stored in the file XDATA.xml.

At that time the problem solved is a 1D rod subjected to a shock. 
The rod has a purely fluid behavior. No elasticity nor plasticity is taken into account.

Two rupture models are available :

1) *imposed pressure* : once the cell reach rupture criterion, the pressure is imposed at the value specified in XDATA.xml)

2) *enrichment* : once the cell reach rupture criterion, kinematic and thermodynamic fields are enriched with an heaviside shape function.

# Future work

* Implement elasticity and plasticity
* Use a non condensed mass matrix to improve results
* Check CFL stability with respect to mass matrix
* Allow multiple discontinuity / enrichment

# Miscellaneous
To launch a profiling under ipython type:
%run -p -r -s module -s cumulative VNR_1D_EnrichedElement.py
