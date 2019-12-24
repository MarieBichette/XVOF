# XVOF
One dimensional hydro code for testing xfem method. 
This software is an experimental work trying to modelize spall rupture thanks to enrichment and partition of unity 
method in a finite volume code. Hansbo & Hansbo enrichment is used.

Input data are stored in the file XDATA.xml.

At that time the problem solved is a 1D rod subjected to a shock.
Material behavior includes simple models for hydrodynamics, elasticity, plasticity.


Two rupture models are available :

1) *imposed pressure* : once the cell reach rupture criterion, the pressure is imposed at the value specified in XDATA.xml

2) *enrichment* : once the cell reach rupture criterion, kinematic and thermodynamic fields are enriched with an heaviside shape function.

Mass matrix can be lumped to improve calculation time.

# Future work
* Check CFL stability with respect to mass matrix

# Miscellaneous
To launch a profiling under ipython type:
%run -p -r -s module -s cumulative VNR_1D_EnrichedElement.py
