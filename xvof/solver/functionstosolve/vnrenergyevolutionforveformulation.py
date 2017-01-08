# -*- coding: iso-8859-1 -*-
"""
Implementing VnrEnergyEvolutionForVolumeEnergyFormulation class
"""

import numpy as np
from xvof.solver.functionstosolve.functiontosolvebase import FunctionToSolveBase


class VnrEnergyEvolutionForVolumeEnergyFormulation(FunctionToSolveBase):
    """
    Defines the evolution function of the internal energy that must vanish in VNR scheme
    [v, e] formulation
    """
    def __init__(self):
        super(VnrEnergyEvolutionForVolumeEnergyFormulation, self).__init__()

    def computeFunctionAndDerivative(self, var_value, mask):
        nrj = var_value[mask]
        eos = self._variables['EquationOfState']
        old_rho = self._variables['OldDensity'][mask]
        new_rho = self._variables['NewDensity'][mask]
        pressure = self._variables['Pressure'][mask]
        old_nrj = self._variables['OldEnergy'][mask]
        p_i = np.zeros(nrj.shape, dtype=np.float64, order='C')
        dpsurde = np.zeros(nrj.shape, dtype=np.float64, order='C')
        dummy = np.zeros(nrj.shape, dtype=np.float64, order='C')
        eos.solveVolumeEnergy(1. / new_rho, nrj, p_i, dummy, dpsurde)
        # Function to vanish
        delta_v = 1. / new_rho - 1. / old_rho
        func = nrj + p_i * delta_v / 2. + pressure * delta_v / 2. - old_nrj
        # Derivative of the function with respect to internal energy
        dfunc = 1 + dpsurde * delta_v / 2.
        return (func, dfunc)
