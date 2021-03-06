# -*- coding: iso-8859-1 -*-
"""
Implementing VnrEnergyEvolutionForVolumeEnergyFormulation class
"""

import numpy as np
from xfv.src.solver.functionstosolve.functiontosolvebase import FunctionToSolveBase


class VnrEnergyEvolutionForVolumeEnergyFormulation(FunctionToSolveBase):
    """
    Defines the evolution function of the internal energy that must vanish in VNR scheme
    [v, e] formulation
    """
    def computeFunctionAndDerivative(self, var_value, mask):
        _mask = np.where(mask)
        nrj = var_value[_mask]
        new_spec_vol = self._variables['NewSpecificVolume'][_mask]

        eos = self._variables['EquationOfState']
        p_i = np.ndarray(nrj.shape, dtype=np.float64, order='C')
        dpsurde = np.ndarray(nrj.shape, dtype=np.float64, order='C')

        eos.solve_volume_energy(new_spec_vol, nrj, p_i, dpsurde)
        # Function to vanish
        delta_v = new_spec_vol - self._variables['OldSpecificVolume'][_mask]
        func = (nrj + (p_i + self._variables['Pressure'][_mask]) * delta_v * 0.5 -
                self._variables['OldEnergy'][_mask])
        # Derivative of the function with respect to internal energy
        dfunc = 1 + dpsurde * delta_v * 0.5
        return func, dfunc
