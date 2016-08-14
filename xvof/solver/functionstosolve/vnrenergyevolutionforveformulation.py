# -*- coding: iso-8859-1 -*-
"""
Implementing VnrEnergyEvolutionForVolumeEnergyFormulation class
"""

from xvof.solver.functionstosolve.functiontosolvebase import FunctionToSolveBase


class VnrEnergyEvolutionForVolumeEnergyFormulation(FunctionToSolveBase):
    """
    Defines the evolution function of the internal energy that must vanish in VNR scheme
    [v, e] formulation
    """
    def __init__(self):
        super(VnrEnergyEvolutionForVolumeEnergyFormulation, self).__init__()

    def computeFunctionAndDerivative(self, newton_variable_value):
        nrj = newton_variable_value
        eos = self._variables['EquationOfState']
        old_rho = self._variables['OldDensity']
        new_rho = self._variables['NewDensity']
        pressure = self._variables['Pressure']
        old_nrj = self._variables['OldEnergy']
        p_i = self._variables['NewPressure']
        dpsurde = self._variables['NewDpOverDe']
        dummy = self._variables['NewSoundSpeed']
        eos.solveVolumeEnergy(1. / new_rho, nrj, p_i, dummy, dpsurde)
        # Function to vanish
        delta_v = 1. / new_rho - 1. / old_rho
        func = nrj + p_i * delta_v / 2. + pressure * delta_v / 2. - old_nrj
        # Derivative of the function with respect to internal energy
        dfunc = 1 + dpsurde * delta_v / 2.
        return (func, dfunc)
