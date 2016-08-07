# -*- coding: iso-8859-1 -*-
"""
Classe définissant un élément en 1d
"""
import ctypes
import numpy as np
import os

from xvof.data.data_container import DataContainer
from xvof.element import Element
from xvof.solver.functionstosolve.vnrenergyevolutionforveformulation import VnrEnergyEvolutionForVolumeEnergyFormulation
from xvof.solver.newtonraphson import NewtonRaphson


class OneDimensionElement(Element):
    """
    Une classe pour les ï¿½lï¿½ments en 1D
    """
    nbr_noeuds = 2

    @classmethod
    def computePseudo(cls, delta_t, rho_old, rho_new, size_new,
                      cel_son, a_pseudo, b_pseudo):
        """
        Calcul de la pseudo
        """
        # pylint: disable=too-many-arguments
        # 8 arguments semblent nï¿½cessaires...
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
        return pseudo

    @classmethod
    def computeTimeStep(cls, cfl, rho_old, rho_new, taille_new, cson_new,
                        pseudo):
        """
        Calcul du pas de temps
        """
        # pylint: disable=too-many-arguments
        # 7 arguments pour cette mï¿½thode cela semble ok
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
        if DataContainer().hasExternalSolver():
            self._external_library = DataContainer().getExternalSolverPath()
        else:
            self._external_library = None
        if self._external_library is not None :
            _path = os.path.join(*(os.path.split(__file__)[:-1] + (self._external_library,)))
            self._mod = ctypes.cdll.LoadLibrary(_path)
            self._computePressureExternal = self._mod.launch_vnr_resolution
            self._computePressureExternal.argtypes = ([ctypes.POINTER(ctypes.c_double), ] * 4 +
                [ctypes.c_int, ] + [ctypes.POINTER(ctypes.c_double), ] * 3)

    # --------------------------------------------------------
    #            DEFINITION DES PROPRIETES                   #
    # --------------------------------------------------------
    @property
    def mass(self):
        """ Masse de l'ï¿½lï¿½ment """
        return self.size_t * self.proprietes.geometric.section * \
            self.density.current_value

    @property
    def pressure_field(self):
        """
        Pressure field
        """
        return self.pressure.current_value

    @property
    def density_field(self):
        """
        Density field
        """
        return self.density.current_value

    @property
    def energy_field(self):
        """
        Internal energy field
        """
        return self.energy.current_value

    @property
    def pseudoviscosity_field(self):
        """
        Pseudoviscosity field
        """
        return self.pseudo.current_value
    # --------------------------------------------------------
    #            DEFINITION DES METHODES                     #
    # --------------------------------------------------------
    def computeNewPressure(self, mask):
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
#         try:
        if self._external_library is not None:
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
            self._computePressureExternal(old_density, new_density, pressure, old_energy, pb_size,
                                          solution, new_pressure, new_vson)
            self.energy.new_value[mask] = solution[0:nbr_cells_to_solve]
            self.pressure.new_value[mask] = new_pressure[0:nbr_cells_to_solve]
            self.sound_velocity.new_value[mask] = new_vson[0:nbr_cells_to_solve]
        else:
            my_variables = {'EquationOfState': self.proprietes.material.eos,
                            'OldDensity': density_current_value,
                            'NewDensity': density_new_value,
                            'Pressure': pressure_current_value + 2. * pseudo_current_value,
                            'OldEnergy': energy_current_value}
            self._function_to_vanish.setVariables(my_variables)
            solution = self._solver.computeSolution(energy_current_value)
            new_pressure_value, _, new_vson_value = \
                self.proprietes.material.eos.solveVolumeEnergy(1. / density_new_value, solution)
            self.energy.new_value[mask] = solution
            self.pressure.new_value[mask] = new_pressure_value
            self.sound_velocity.new_value[mask] = new_vson_value
            self._function_to_vanish.eraseVariables()
#         except ValueError as err:
#             raise err

    def computeSize(self, topologie, vecteur_coord_noeuds):
        """
        Calcul de la longueur de l'ï¿½lï¿½ment (ï¿½ t)
        """
        for ielem in xrange(self._shape[0]):
            ind_nodes = topologie.getNodesBelongingToCell(ielem)
            self._size_t[ielem] = abs(vecteur_coord_noeuds[ind_nodes[0]] -
                                      vecteur_coord_noeuds[ind_nodes[1]])

    def computeNewSize(self, topologie, vecteur_coord_noeuds, mask, time_step=None):
        """
        Calcul de la nouvelle longueur de l'élément (à t+dt)
        """
        connectivity = np.array(topologie._nodes_belonging_to_cell)
        size_t_plus_dt = abs(vecteur_coord_noeuds[connectivity[:, 0]] -
                             vecteur_coord_noeuds[connectivity[:, 1]]).reshape(self.number_of_cells)
        self._size_t_plus_dt[mask] = size_t_plus_dt[mask]

    def computeNewDensity(self, mask):
        """
        Calcul de la densite a l'instant t+dt basee sur
        la conservation de la masse
        """
        self.density.new_value[mask] = \
            self.density.current_value[mask] * self.size_t[mask] / self.size_t_plus_dt[mask]

    def computeNewPseudo(self, delta_t, mask):
        """
        Calcul de la nouvelle pseudo
        """
        self.pseudo.new_value[mask] = \
            OneDimensionElement.computePseudo(delta_t, self.density.current_value[mask], self.density.new_value[mask],
                                              self.size_t_plus_dt[mask], self.sound_velocity.current_value[mask],
                                              self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

    def computeNewTimeStep(self, mask):
        """
        Calcul du pas de temps dans l'ï¿½lï¿½ment
        """
        cfl = self.proprietes.numeric.cfl
        dt = \
            OneDimensionElement.computeTimeStep(cfl, self.density.current_value, self.density.new_value,
                                                self.size_t_plus_dt, self.sound_velocity.new_value,
                                                self.pseudo.current_value)
        self._dt[mask] = dt[mask]

    def imposePressure(self, ind_cell, pression):
        """
        On impose la pression ï¿½ t+dt (par exemple pour endommagement)
        """
        self.pressure.new_value[ind_cell] = pression
