# -*- coding: iso-8859-1 -*-
"""
Implementation of the OneDimensionCell class
"""
import ctypes
import numpy as np
import os

from xvof.cell import Cell
from xvof.data.data_container import DataContainer
from xvof.solver.functionstosolve.vnrenergyevolutionforveformulation import VnrEnergyEvolutionForVolumeEnergyFormulation
from xvof.solver.newtonraphson import NewtonRaphson


class OneDimensionCell(Cell):
    """
    A class for one dimension cells
    """
    @classmethod
    def compute_pseudo(cls, delta_t, rho_old, rho_new, size_new, cel_son, a_pseudo, b_pseudo):
        """
        Computation of artificial viscosity
        """
        # pylint: disable=too-many-arguments
        # 8 arguments semblent nécessaires
        vnt = 1. / rho_old
        vnplusun = 1. / rho_new
        vnplusundemi = 0.5 * (vnt + vnplusun)
        vpointnplusundemi = 1. / delta_t * (vnplusun - vnt)
        divu = vpointnplusundemi / vnplusundemi
        pseudo = np.zeros(rho_old.shape, dtype=np.float64, order='C')
        mask = divu < 0.
        pseudo[mask] = (1. / vnplusundemi[mask] *
                        (a_pseudo * size_new[mask] ** 2 * vpointnplusundemi[mask] ** 2 / vnplusundemi[mask] ** 2 +
                         b_pseudo * size_new[mask] * cel_son[mask] * abs(vpointnplusundemi[mask]) / vnplusundemi[mask]))
        return pseudo

    @classmethod
    def compute_time_step(cls, cfl, cfl_pseudo, rho_old, rho_new, taille_new, cson_new, pseudo_old, pseudo_new):
        """
        Computation of the time step
        """
        # pylint: disable=too-many-arguments
        # 7 arguments pour cette mï¿½thode cela semble ok
        local_cson = np.copy(cson_new) ** 2
        mask_q = pseudo_new != 0.
        drho = np.abs((rho_new - rho_old) / rho_old)
        dpseudo = (pseudo_new - pseudo_old)
        mask_r = drho > 1.e-04
        mask_local_cson = np.logical_and(mask_q, mask_r)
        pseudo_sound_speep_square = np.abs(cfl_pseudo * dpseudo[mask_local_cson] /
                                           (rho_new[mask_local_cson] - rho_old[mask_local_cson]))
        local_cson[mask_local_cson] += pseudo_sound_speep_square
        local_cson **= 0.5
        delta_t = cfl * taille_new / local_cson
        return delta_t

    def __init__(self, number_of_elements):
        Cell.__init__(self, number_of_elements)
        self._function_to_vanish = VnrEnergyEvolutionForVolumeEnergyFormulation()
        self._solver = NewtonRaphson(self._function_to_vanish)
        self._data_path_file = os.path.join(os.path.curdir, "XDATA.xml")
        #
        if DataContainer(self._data_path_file).hasExternalSolver():
            self._external_library = DataContainer(self._data_path_file).getExternalSolverPath()
        else:
            self._external_library = None
        if self._external_library is not None:
            _path = os.path.join(*(os.path.split(__file__)[:-1] + (self._external_library,)))
            self._mod = ctypes.cdll.LoadLibrary(_path)
            self._computePressureExternal = self._mod.launch_vnr_resolution
            self._computePressureExternal.argtypes = ([ctypes.POINTER(ctypes.c_double), ] * 4 +
                                                      [ctypes.c_int, ] + [ctypes.POINTER(ctypes.c_double), ] * 3)

    def compute_mass(self):
        """
        Compute mass of the cells
        """
        self._mass = self.size_t * DataContainer(self._data_path_file).geometric.section * self.density.current_value

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
    def artificial_viscosity_field(self):
        """
        Pseudoviscosity field
        """
        return self.pseudo.current_value

    def _compute_new_pressure_with_external_lib(self, density_current, density_new, pressure_current, pseudo_current,
                                                energy_current, energy_new, pressure_new, vson_new):
        """
        Computation of the set (internal energy, pressure, sound velocity) for v-e formulation thanks
        to external C library
        """
        pb_size = ctypes.c_int()
        pb_size.value = energy_new.shape[0]
        c_density = density_current.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        n_density = density_new.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        true_pressure = (pressure_current + 2. * pseudo_current)
        c_pressure = true_pressure.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_energy = energy_current.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        n_energy = energy_new.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        n_pressure = pressure_new.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        n_sound_speed = vson_new.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        self._computePressureExternal(c_density, n_density, c_pressure, c_energy, pb_size,
                                      n_energy, n_pressure, n_sound_speed)
        energy_n = n_energy[0:pb_size.value]
        pressure_n = n_pressure[0:pb_size.value]
        vson_n = n_sound_speed[0:pb_size.value]
        return energy_n, pressure_n, vson_n

    def compute_new_pressure(self, mask):
        """
        Computation of the set (internal energy, pressure, sound velocity) for v-e formulation
        """
        if self._external_library is not None:
            self.energy.new_value[mask], self.pressure.new_value[mask], self.sound_velocity.new_value[mask] = \
                self._compute_new_pressure_with_external_lib(self.density.current_value[mask],
                                                             self.density.new_value[mask],
                                                             self.pressure.current_value[mask],
                                                             self.pseudo.current_value[mask],
                                                             self.energy.current_value[mask],
                                                             self.energy.new_value[mask],
                                                             self.pressure.new_value[mask],
                                                             self.sound_velocity.new_value[mask])
        else:
            my_variables = {'EquationOfState': DataContainer(self._data_path_file).material.eos,
                            'OldDensity': self.density.current_value[mask],
                            'NewDensity': self.density.new_value[mask],
                            'Pressure': self.pressure.current_value[mask] + 2. * self.pseudo.current_value[mask],
                            'OldEnergy': self.energy.current_value[mask]}
            self._function_to_vanish.setVariables(my_variables)
            self.energy.new_value[mask] = self._solver.computeSolution(self.energy.current_value[mask])
            # Eos call to determine final pressure and sound speed values
            shape = self.energy.new_value[mask].shape
            new_pressure_value = np.zeros(shape, dtype=np.float64, order='C')
            new_vson_value = np.zeros(shape, dtype=np.float64, order='C')
            dummy = np.zeros(shape, dtype=np.float64, order='C')
            my_variables['EquationOfState'].solveVolumeEnergy(
                    1./ my_variables['NewDensity'], self.energy.new_value[mask], new_pressure_value, new_vson_value,
                    dummy)
            self.pressure.new_value[mask] = new_pressure_value
            self.sound_velocity.new_value[mask] = new_vson_value
            self._function_to_vanish.eraseVariables()

    def compute_size(self, topologie, vecteur_coord_noeuds):
        """
        Computation of the cells initial length
        """
        connectivity = topologie.nodes_belonging_to_cell
        self._size_t = abs(vecteur_coord_noeuds[connectivity[:,0]] - vecteur_coord_noeuds[connectivity[:,1]]).flatten()

    def compute_new_size(self, topologie, vecteur_coord_noeuds, mask):
        """
        Computation of the cells length at time t+dt
        """
        connectivity = topologie.nodes_belonging_to_cell
        size_t_plus_dt = abs(vecteur_coord_noeuds[connectivity[:,0]] - vecteur_coord_noeuds[connectivity[:,1]])
        self._size_t_plus_dt[mask] = size_t_plus_dt[mask].flatten()

    def compute_new_density(self, mask):
        """
        Computation of the density of the cells at time t+dt using mass conservation principle
        """
        self.density.new_value[mask] =  self.density.current_value[mask] * self.size_t[mask] / self.size_t_plus_dt[mask]

    def compute_new_pseudo(self, delta_t, mask):
        """
        Computation of cells artificial viscosity at time t+dt/2
        """
        self.pseudo.new_value[mask] = OneDimensionCell.compute_pseudo(delta_t, self.density.current_value[mask],
                                                                      self.density.new_value[mask],
                                                                      self.size_t_plus_dt[mask],
                                                                      self.sound_velocity.current_value[mask],
                                                                      DataContainer(self._data_path_file).numeric.a_pseudo,
                                                                      DataContainer(self._data_path_file).numeric.b_pseudo)

    def compute_new_time_step(self, mask):
        """
        Computation of the time step in the cells at time t+dt
        """
        cfl = DataContainer(self._data_path_file).numeric.cfl
        cfl_pseudo = DataContainer(self._data_path_file).numeric.cfl_pseudo
        dt = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, self.density.current_value, self.density.new_value,
                                                self.size_t_plus_dt, self.sound_velocity.new_value,
                                                self.pseudo.current_value, self.pseudo.new_value)
        self._dt[mask] = dt[mask]

    def impose_pressure(self, ind_cell, pression):
        """
        Pressure imposition
        """
        self.pressure.new_value[ind_cell] = pression
