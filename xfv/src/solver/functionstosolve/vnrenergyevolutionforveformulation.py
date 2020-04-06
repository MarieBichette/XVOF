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
    def computeFunctionAndDerivative(self, var_value, mask):  # pylint: disable=too-many-locals
        _mask = np.where(mask)
        nrj = var_value[_mask]
        eos = self._variables['EquationOfState']
        old_rho = self._variables['OldDensity'][_mask]
        new_rho = self._variables['NewDensity'][_mask]
        new_spec_vol = 1. / new_rho
        pressure = self._variables['Pressure'][_mask]
        old_nrj = self._variables['OldEnergy'][_mask]
        p_i = np.ndarray(nrj.shape, dtype=np.float64, order='C')
        dpsurde = np.ndarray(nrj.shape, dtype=np.float64, order='C')
        eos.solve_volume_energy(new_spec_vol, nrj, p_i, dpsurde)
        # Function to vanish
        delta_v = new_spec_vol - 1. / old_rho
        func = nrj + p_i * delta_v / 2. + pressure * delta_v / 2. - old_nrj
        # Derivative of the function with respect to internal energy
        dfunc = 1 + dpsurde * delta_v / 2.
        return func, dfunc
