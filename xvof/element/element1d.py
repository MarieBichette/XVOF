#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe d�finissant un �l�ment en 1d
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
    Une classe pour les �l�ments en 1D
    """
    nbr_noeuds = 2

    @classmethod
    def computePseudo(cls, delta_t, rho_old, rho_new, size_new,
                      cel_son, a_pseudo, b_pseudo):
        """
        Calcul de la pseudo
        """
        # pylint: disable=too-many-arguments
        # 8 arguments semblent n�cessaires...
        vnt = 1. / rho_old
        vnplusun = 1. / rho_new
        vnplusundemi = 0.5 * (vnt + vnplusun)
        vpointnplusundemi = 1. / delta_t * (vnplusun - vnt)
        divu = vpointnplusundemi / vnplusundemi
        pseudo = np.zeros(rho_old.shape, dtype=np.float64, order='C')
        mask = divu < 0.
        pseudo[mask] = 1. / vnplusundemi[mask] * \
            (a_pseudo * size_new[mask] ** 2 * vpointnplusundemi[mask] ** 2 / vnplusundemi[mask] ** 2 +
             b_pseudo * size_new[mask] * cel_son[mask] *
             abs(vpointnplusundemi[mask]) / vnplusundemi[mask])
#         for icell in xrange(rho_old.shape[0]):
#             if divu[icell] < 0.:
#                 pseudo[icell] = 1. / vnplusundemi[icell] * \
#                     (a_pseudo * size_new[icell] ** 2 * vpointnplusundemi[icell] ** 2 / vnplusundemi[icell] ** 2 +
#                      b_pseudo * size_new[icell] * cel_son[icell] *
#                      abs(vpointnplusundemi[icell]) / vnplusundemi[icell])
        return pseudo

    @classmethod
    def computeTimeStep(cls, cfl, rho_old, rho_new, taille_new, cson_new,
                        pseudo):
        """
        Calcul du pas de temps
        """
        # pylint: disable=too-many-arguments
        # 7 arguments pour cette m�thode cela semble ok
        delta_t = np.zeros(rho_old.shape, dtype=np.float64, order='C')
        drho = rho_new - rho_old
        mask = drho > 0.1
        delta_t[mask] = cfl * taille_new[mask] / (cson_new[mask] ** 2 + 2. * pseudo[mask] / drho[mask]) ** 0.5
        mask = drho <= 0.1
        delta_t[mask] = cfl * taille_new[mask] / cson_new[mask]
        return delta_t

    def __init__(self, number_of_elements, proprietes):
        Element.__init__(self, number_of_elements, proprietes)
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
        """ Masse de l'�l�ment """
        return self.size_t * self.proprietes.geometric.section * \
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
        density_current_value = self.density.current_value[mask]
        density_new_value = self.density.new_value[mask]
        pressure_current_value = self.pressure.current_value[mask]
        pseudo_current_value = self.pseudo.current_value[mask]
        energy_current_value = self.energy.current_value[mask]
        energy_new_value = self.energy.new_value[mask]
        shape = energy_new_value.shape
        nbr_cells_to_solve = shape[0]
        solution_value = np.zeros(shape, dtype=np.float64, order='C')
        new_pressure_value = np.zeros(shape, dtype=np.float64, order='C')
        new_vson_value = np.zeros(shape, dtype=np.float64, order='C')
        try:
            if EXTERNAL_LIBRARY is not None:
                pb_size = ctypes.c_int()
                #
                old_density = density_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                new_density = density_new_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                tmp = (pressure_current_value + 2. * pseudo_current_value)
                pressure = tmp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                old_energy = energy_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                pb_size.value = nbr_cells_to_solve
                solution = solution_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                new_pressure = new_pressure_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                new_vson = new_vson_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
#                 if nbr_cells_to_solve != 1001:
#                     import ipdb; ipdb.set_trace()
                self._computePressureExternal(old_density, new_density, pressure, old_energy, pb_size,
                                              solution, new_pressure, new_vson)
                self.energy.new_value[mask] = solution[0:nbr_cells_to_solve]
                self.pressure.new_value[mask] = new_pressure[0:nbr_cells_to_solve]
                self.sound_velocity.new_value[mask] = new_vson[0:nbr_cells_to_solve]
#                print "{:15.9g} | {:15.9g} | {:15.9g}".format(self.energy.new_value, self.pressure.new_value, self.sound_velocity.new_value)
            else:
                my_variables = {'EquationOfState': self.proprietes.material.eos,
                                'OldDensity': density_current_value,
                                'NewDensity': density_new_value,
                                'Pressure': pressure_current_value + 2. * pseudo_current_value,
                                'OldEnergy': energy_current_value}
                self._function_to_vanish.setVariables(my_variables)
                self.energy.new_value = self._solver.computeSolution(energy_current_value)
                for icell in xrange(self.number_of_cells):
                    self.pressure.new_value[icell], _, self.sound_velocity.new_value[icell] = \
                        self.proprietes.material.eos.solveVolumeEnergy(1. / density_new_value[icell], energy_new_value[icell])
                self._function_to_vanish.eraseVariables()
#            print "{:15.9g} | {:15.9g} | {:15.9g}".format(self.energy.new_value, self.pressure.new_value, self.sound_velocity.new_value)
#            raw_input()
        except ValueError as err:
            raise err

    def computeSize(self, topologie, vecteur_coord_noeuds):
        """
        Calcul de la longueur de l'�l�ment (� t)
        """
        for ielem in xrange(self._shape[0]):
            ind_nodes = topologie.getNodesBelongingToCell(ielem)
            self._size_t[ielem] = abs(vecteur_coord_noeuds[ind_nodes[0]] -
                                      vecteur_coord_noeuds[ind_nodes[1]])

    def computeNewSize(self, topologie, vecteur_coord_noeuds, time_step=None):
        """
        Calcul de la nouvelle longueur de l'�l�ment (� t+dt)
        """
        connectivity = np.array(topologie._nodes_belonging_to_cell)
        self._size_t_plus_dt = abs(vecteur_coord_noeuds[connectivity[:, 0]] -
                                   vecteur_coord_noeuds[connectivity[:, 1]]).reshape(self.number_of_cells)

    def computeNewDensity(self):
        """
        Calcul de la densit� � l'instant t+dt bas� sur
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
                                    self.size_t_plus_dt, self.sound_velocity.current_value,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

    def computeNewTimeStep(self):
        """
        Calcul du pas de temps dans l'�l�ment
        """
        cfl = self.proprietes.numeric.cfl
        self._dt = \
            Element1d.computeTimeStep(cfl, self.density.current_value, self.density.new_value,
                                      self.size_t_plus_dt, self.sound_velocity.new_value,
                                      self.pseudo.current_value)

    def imposePressure(self, ind_cell, pression):
        """
        On impose la pression � t+dt (par exemple pour endommagement)
        """
        self.pressure.new_value[ind_cell] = pression
