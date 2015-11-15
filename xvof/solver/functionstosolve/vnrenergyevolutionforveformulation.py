#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant la fonction d'évolution de l'énergie interne à annuler dans le schéma VNR
Formulation V-E
"""
import numpy as np
from xvof.solver.functionstosolve.functiontosolvebase import FunctionToSolveBase


class VnrEnergyEvolutionForVolumeEnergyFormulation(FunctionToSolveBase):
    '''
    Classe définissant la fonction d'évolution de l'énergie interne à annuler dans le schéma VNR
    Formulation V-E
    '''
    def __init__(self):
        super(VnrEnergyEvolutionForVolumeEnergyFormulation, self).__init__()

    def computeFunctionAndDerivative(self, newton_variable_value):
        nrj = newton_variable_value
        eos = self._variables['EquationOfState']
        old_rho = self._variables['OldDensity']
        new_rho = self._variables['NewDensity']
        pressure = self._variables['Pressure']
        old_nrj = self._variables['OldEnergy']
        p_i = np.zeros(old_rho.shape, dtype=np.float64, order='C')
        dpsurde = np.zeros(old_rho.shape, dtype=np.float64, order='C')
        for icell in xrange(old_rho.shape[0]):
            (p_i[icell], dpsurde[icell], dummy) = eos.solveVolumeEnergy(1. / new_rho[icell], nrj[icell])
        # Fonction à annuler
        delta_v = 1. / new_rho - 1. / old_rho
        func = nrj + p_i * delta_v / 2. + pressure * delta_v / 2. - old_nrj
        # Dérivée de la fonction à annuler
        dfunc = 1 + dpsurde * delta_v / 2.
        return (func, dfunc)
