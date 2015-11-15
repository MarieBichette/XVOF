#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un élément en 1d
"""
# --------------------------------------------------------
#               IMPORTATIONS DIVERSES                    #
# --------------------------------------------------------
import ctypes
import os
import numpy as np
import numpy.ma as ma

from xvof.element import Element
from xvof.solver.functionstosolve.vnrenergyevolutionforveformulation import VnrEnergyEvolutionForVolumeEnergyFormulation
from xvof.solver.newtonraphson import NewtonRaphson


EXTERNAL_LIBRARY = 'vnr_internal_energy_evolution.so'
# EXTERNAL_LIBRARY = None
# --------------------------------------------------------
#        DEFINITION DES CLASSES ET FONCTIONS             #
# --------------------------------------------------------
class Element1d(Element):
    """
    Une classe pour les éléments en 1D
    """
    nbr_noeuds = 2

    @classmethod
    def computePseudo(cls, delta_t, rho_old, rho_new, size_new,
                      cel_son, a_pseudo, b_pseudo):
        """
        Calcul de la pseudo
        """
        # pylint: disable=too-many-arguments
        # 8 arguments semblent nécessaires...
        vnt = 1. / rho_old
        vnplusun = 1. / rho_new
        vnplusundemi = 0.5 * (vnt + vnplusun)
        vpointnplusundemi = 1. / delta_t * (vnplusun - vnt)
        divu = vpointnplusundemi / vnplusundemi
        pseudo = 0.
        if divu < 0.:
            pseudo = 1. / vnplusundemi * \
                (a_pseudo * size_new ** 2 * vpointnplusundemi ** 2 / vnplusundemi ** 2 +
                 b_pseudo * size_new * cel_son *
                 abs(vpointnplusundemi) / vnplusundemi)
        return pseudo

    @classmethod
    def computeTimeStep(cls, cfl, rho_old, rho_new, taille_new, cson_new,
                        pseudo):
        """
        Calcul du pas de temps
        """
        # pylint: disable=too-many-arguments
        # 7 arguments pour cette méthode cela semble ok
        delta_t = 0.
        if (rho_new - rho_old) > 0.1:
            delta_t = cfl * taille_new / ((cson_new ** 2 + 2. * pseudo /
                                           (rho_new - rho_old)) ** 0.5)
        else:
            delta_t = cfl * taille_new / cson_new
        return delta_t

    def __init__(self, number_of_elements, proprietes):
        Element.__init__(self, number_of_elements, proprietes, pressure_offset=2)
        self._function_to_vanish = VnrEnergyEvolutionForVolumeEnergyFormulation()
        self._solver = NewtonRaphson(self._function_to_vanish)
        #
        if EXTERNAL_LIBRARY is not None :
            _path = os.path.join(*(os.path.split(__file__)[:-1] + (EXTERNAL_LIBRARY,)))
            self._mod = ctypes.cdll.LoadLibrary(_path)
            self._computePressureExternal = self._mod.launch_vnr_resolution
            self._computePressureExternal.argtypes = ([ctypes.POINTER(ctypes.c_double), ] * 4 +
                [ctypes.c_int, ] + [ctypes.POINTER(ctypes.c_double), ] * 3)

    # --------------------------------------------------------
    #            DEFINITION DES PROPRIETES                   #
    # --------------------------------------------------------
    @property
    def masse(self):
        """ Masse de l'élément """
        return self.taille_t * self.proprietes.geometric.section * \
            self.density.current_value

    # --------------------------------------------------------
    #            DEFINITION DES METHODES                     #
    # --------------------------------------------------------
    def computeNewPressure(self, mask=None):
        """
        Calcul du triplet energie, pression, vitesse du son
        au pas de temps suivant
        Formulation v-e
        """
        density_current_value = ma.masked_array(self.density.current_value)
        density_new_value = ma.masked_array(self.density.new_value)
        pressure_current_value = ma.masked_array(self.pressure.current_value)
        pseudo_current_value = ma.masked_array(self.pseudo.current_value)
        energy_current_value = ma.masked_array(self.energy.current_value)
        energy_new_value = ma.masked_array(self.energy.new_value)
        if mask is not None:
            density_current_value[mask] = ma.masked
            density_new_value[mask] = ma.masked
            pressure_current_value[mask] = ma.masked
            pseudo_current_value[mask] = ma.masked
            energy_current_value[mask] = ma.masked
        try:
            if EXTERNAL_LIBRARY is not None:
                old_density = ctypes.c_double()
                new_density = ctypes.c_double()
                pressure = ctypes.c_double()
                old_energy = ctypes.c_double()
                pb_size = ctypes.c_int()
                solution = ctypes.c_double()
                new_pressure = ctypes.c_double()
                new_vson = ctypes.c_double()
                #
                old_density.value = density_current_value
                new_density.value = density_new_value
                pressure.value = pressure_current_value + 2. * pseudo_current_value
                old_energy.value = energy_current_value
                pb_size.value = 1
                self._computePressureExternal(old_density, new_density, pressure, old_energy, pb_size, solution,
                                              new_pressure, new_vson)
                self.energy.new_value = solution.value
                self.pressure.new_value = new_pressure.value
                self.sound_velocity.new_value = new_vson.value
#                print "{:15.9g} | {:15.9g} | {:15.9g}".format(self.energy.new_value, self.pressure.new_value, self.sound_velocity.new_value)
            else:
                my_variables = {'EquationOfState': self.proprietes.material.eos,
                                'OldDensity': density_current_value,
                                'NewDensity': density_new_value,
                                'Pressure': pressure_current_value + 2. * pseudo_current_value,
                                'OldEnergy': energy_current_value}
                self._function_to_vanish.setVariables(my_variables)
                self.energy.new_value = self._solver.computeSolution(energy_current_value)
                self.pressure.new_value, _, self.sound_velocity.new_value = \
                    self.proprietes.material.eos.solveVolumeEnergy(1. / density_new_value, energy_new_value)
                self._function_to_vanish.eraseVariables()
#            print "{:15.9g} | {:15.9g} | {:15.9g}".format(self.energy.new_value, self.pressure.new_value, self.sound_velocity.new_value)
#            raw_input()
        except ValueError as err:
            print "Element concerné : {}".format(self)
            raise err

    def computeSize(self, topologie, vecteur_coord_noeuds):
        """
        Calcul de la longueur de l'élément (à t)
        """
        for ielem in xrange(self._shape[0]):
            ind_nodes = topologie.getNodesBelongingToCell(ielem)
            self._size_t[ielem] = abs(vecteur_coord_noeuds[ind_nodes[0]] -
                                      vecteur_coord_noeuds[ind_nodes[1]])

    def computeNewSize(self, topologie, vecteur_coord_noeuds, time_step=None):
        """
        Calcul de la nouvelle longueur de l'élément (à t+dt)
        """
        for ielem in xrange(self._shape[0]):
            ind_nodes = topologie.getNodesBelongingToCell(ielem)
            self._size_t_plus_dt[ielem] = abs(vecteur_coord_noeuds[ind_nodes[0]] -
                                              vecteur_coord_noeuds[ind_nodes[1]])

    def computeNewDensity(self):
        """
        Calcul de la densité à l'instant t+dt basé sur
        la conservation de la masse
        """
        self.density.new_value = \
            self.density.current_value * self.size_t / self.size_t_plus_dt

    def computeNewPseudo(self, delta_t):
        """
        Calcul de la nouvelle pseudo
        """
        self.pseudo.new_value = \
            Element1d.computePseudo(delta_t, self.density.current_value, self.density.new_value,
                                    self.taille_t_plus_dt, self.sound_velocity.current_value,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

    def computeNewTimeStep(self):
        """
        Calcul du pas de temps dans l'élément
        """
        cfl = self.proprietes.numeric.cfl
        self._dt = \
            Element1d.computeTimeStep(cfl, self.density.current_value, self.density.new_value,
                                      self.taille_t_plus_dt, self.sound_velocity.new_value,
                                      self.pseudo.current_value)

    def imposePressure(self, ind_cell, pression):
        """
        On impose la pression à t+dt (par exemple pour endommagement)
        """
        self.pressure.new_value[ind_cell] = pression
