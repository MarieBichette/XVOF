#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe d�finissant la fonction d'�volution de l'�nergie interne � annuler dans le sch�ma VNR
Formulation V-E
"""
from xvof.solver.functionstosolve.functiontosolvebase import FunctionToSolveBase


class VnrEnergyEvolutionForVolumeEnergyFormulation(FunctionToSolveBase):
    '''
    Classe d�finissant la fonction d'�volution de l'�nergie interne � annuler dans le sch�ma VNR
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
        (p_i, dpsurde, dummy) = eos.solveVolumeEnergy(1. / new_rho, nrj)
        # Fonction � annuler
        delta_v = 1. / new_rho - 1. / old_rho
        func = nrj + p_i * delta_v / 2. + pressure * delta_v / 2. - \
            old_nrj
        # D�riv�e de la fonction � annuler
        dfunc = 1 + dpsurde * delta_v / 2.
        return (func, dfunc)