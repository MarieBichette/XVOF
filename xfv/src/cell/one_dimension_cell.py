#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementation of the OneDimensionCell class
"""
import ctypes
import numpy as np
import os

from xfv.src.cell import Cell
from xfv.src.solver.functionstosolve.vnrenergyevolutionforveformulation import VnrEnergyEvolutionForVolumeEnergyFormulation
from xfv.src.solver.newtonraphson import NewtonRaphson
from xfv.src.data.data_container import DataContainer
from xfv.src.utilities.stress_invariants_calculation import compute_J2


# noinspection PyArgumentList
class OneDimensionCell(Cell):
    """
    A class for one dimension cells
    """

    @classmethod
    def apply_equation_of_state(cls, cell, eos, density, density_new, pressure, pressure_new,
                                energy, energy_new, pseudo, cson_new):

        if cell._external_library is not None:
            energy_new_value, pressure_new_value, sound_velocity_new_value = \
                (np.array(x) for x in cell._compute_new_pressure_with_external_lib(
                    density, density_new, pressure, pseudo, energy, energy_new,
                    pressure_new, cson_new))
        else:
            my_variables = {'EquationOfState': eos,
                            'OldDensity': density,
                            'NewDensity': density_new,
                            'Pressure': (pressure + 2. * pseudo),
                            'OldEnergy': energy}
            cell._function_to_vanish.setVariables(my_variables)
            energy_new_value = cell._solver.compute_solution(energy)

            # Eos call to determine final pressure and sound speed values
            shape = energy_new.shape
            pressure_new_value = np.zeros(shape, dtype=np.float64, order='C')
            sound_velocity_new_value = np.zeros(shape, dtype=np.float64, order='C')
            dummy = np.zeros(shape, dtype=np.float64, order='C')
            my_variables['EquationOfState'].solveVolumeEnergy(
                1. / my_variables['NewDensity'], energy_new_value, pressure_new_value,
                sound_velocity_new_value, dummy)
            cell._function_to_vanish.eraseVariables()
        return energy_new_value, pressure_new_value, sound_velocity_new_value

    @classmethod
    def add_elastic_energy_method(cls, dt, density_current, density_new,
                                  stress_dev_current, stress_dev_new, strain_rate_dev):
        """
        Take into account the additional term in internal energy due to elasticity
        """
        # Le 1/2 de la moyenne density_current et density_new se simplifie avec
        # la moyenne stress_dev_current et stress_dev_new
        energy_new_value = dt / (density_new + density_current)\
                           * ((stress_dev_new[:, 0]
                               + stress_dev_current[:, 0]) * strain_rate_dev[:, 0] +
                              (stress_dev_new[:, 1]
                               + stress_dev_current[:, 1]) * strain_rate_dev[:, 1] +
                              (stress_dev_new[:, 2]
                               + stress_dev_current[:, 2]) * strain_rate_dev[:, 2])
        return energy_new_value

    @classmethod
    def general_method_deviator_strain_rate(cls, mask, dt, x_new, u_new):
        """
        Compute the deviateur du tenseur taux de dï¿½formation (calculï¿½ au centre la maille)
        ï¿½ partir des grandeurs vitesse et position au centre des mailles
        :param mask : tableau de booleen pour identifier les cell ï¿½ calculer
        :param dt : time step (float)
<<<<<<< HEAD
        :param x_new : array des coordonnées à l'intant n+1 des neouds à gauche et
        à droite de la cell calculée
        :param u_new : array des vitesses à l'intant n+1 des neouds à gauche et
        à droite de la cell calculée
=======
        :param x_new : array des coordonnï¿½es ï¿½ l'intant n+1 des neouds ï¿½ gauche et ï¿½ droite de la cell calculï¿½e
        :param u_new : array des vitesses ï¿½ l'intant n+1 des neouds ï¿½ gauche et ï¿½ droite de la cell calculï¿½e
>>>>>>> 6b55f92392fa87d4a7349199dbad65dc2e1c3323
        x_new, u_new sont obligatoirement tous de taille (taille de mask, 2)
        """
        strain_rate_dev = np.zeros([u_new.shape[0], 3])
        # Calcul du dï¿½viateur de D
        x_demi = x_new - dt/2. * u_new
        D = (u_new[mask, 1] - u_new[mask, 0]) / (x_demi[mask, 1] - x_demi[mask, 0])  # Dxx
        strain_rate_dev[mask, 0] = 2. / 3. * D
        strain_rate_dev[mask, 1] = - 1. / 3. * D
        strain_rate_dev[mask, 2] = - 1. / 3. * D
        return strain_rate_dev[mask]

    @classmethod
    def compute_pseudo(cls, delta_t, rho_old, rho_new, size_new, cel_son, a_pseudo, b_pseudo):
        """
        Computation of artificial viscosity
        """
        # pylint: disable=too-many-arguments
        # 8 arguments semblent nï¿½cessaires
        vn = 1. / rho_old
        vnplusun = 1. / rho_new
        vnplusundemi = (vn + vnplusun) / 2.
        vpoint = (vnplusun - vn) / delta_t
        divu = vpoint / vnplusundemi
        pseudo = np.zeros(rho_old.shape, dtype=np.float64, order='C')
        mask = np.where(divu < 0.)
        pseudo[mask] = 1. / vnplusundemi[mask] * (
                a_pseudo * size_new[mask] ** 2 * divu[mask] ** 2 +
                b_pseudo * size_new[mask] * cel_son[mask] * abs(divu[mask]))
        return pseudo

    @classmethod
    def compute_time_step(cls, cfl, cfl_pseudo, rho_old, rho_new, taille_new, cson_new,
                          pseudo_old, pseudo_new):
        """
        Computation of the time step
        :param cfl : nombre cfl
        :param cfl_pseudo : cfl par rapport ï¿½ la condition de traitement du choc
        :param rho_old : density at time t
        :param rho_new : density at time t+dt
        :param taille_new : size of element
        :param cson_new : sound velocity at time t+dt
        :param pseudo_old : artificial viscosity at time t
        :param pseudo_new : artificial viscosity at timet+dt
        """
        # pylint: disable=too-many-arguments
        local_cson = np.copy(cson_new) ** 2
        mask_q = pseudo_new != 0.
        drho = np.abs((rho_new - rho_old) / rho_old)
        dpseudo = (pseudo_new - pseudo_old)
        mask_r = drho > 1.e-04
        mask_local_cson = np.logical_and(mask_q, mask_r)
        pseudo_sound_speed_square = np.abs(cfl_pseudo * dpseudo[mask_local_cson] /
                                           (rho_new[mask_local_cson] - rho_old[mask_local_cson]))
        local_cson[mask_local_cson] += pseudo_sound_speed_square
        local_cson **= 0.5
        delta_t = cfl * taille_new / local_cson
        return delta_t

    def __init__(self, number_of_elements):
        super(OneDimensionCell, self).__init__(number_of_elements)

        # By default :all cells are classical (non enriched)
        self._classical = np.ones([number_of_elements, ], dtype=np.bool, order='C')
        self._enrichment_not_concerned = np.ones([number_of_elements, ], dtype=np.bool, order='C')

        self.dtc = []  # pour enregistrement du pas de temps critique

        # elasticity / plasticity:
        self._deviatoric_stress_new = np.zeros([number_of_elements, 3], dtype=np.float64,
                                               order='C')  # dev sxx, syy, szz
        self._deviatoric_stress_current = np.zeros([number_of_elements, 3], dtype=np.float64,
                                                   order='C')  # dev sxx, syy, szz
        self._deviatoric_strain_rate = np.zeros([number_of_elements, 3], dtype=np.float64,
                                                order='C')
        self._equivalent_plastic_strain_rate = np.zeros([number_of_elements, ], dtype=np.float64,
                                                        order='C')
        self._plastic_strain_rate = np.zeros([number_of_elements, 3], dtype=np.float64,
                                             order='C')  # Dp_xx, Dp_yy, Dp_zz

        # Endommagement / CZM
        self._damage_variable = np.zeros([number_of_elements, ], dtype=np.float64, order='C')

        # Solveur pour EOS :
        self._function_to_vanish = VnrEnergyEvolutionForVolumeEnergyFormulation()
        self._solver = NewtonRaphson(self._function_to_vanish)

        if DataContainer().hasExternalSolver():
            self._external_library = DataContainer().getExternalSolverPath()
        else:
            self._external_library = None
        if self._external_library is not None:
            _path = os.path.join(*(os.path.split(__file__)[:-1] + (self._external_library,)))
            self._mod = ctypes.cdll.LoadLibrary(_path)
            self._computePressureExternal = self._mod.launch_vnr_resolution
            self._computePressureExternal.argtypes = ([ctypes.POINTER(ctypes.c_double), ] * 4 +
                                                      [ctypes.c_int, ] +
                                                      [ctypes.POINTER(ctypes.c_double), ] * 3)

    def compute_mass(self):
        """
        Compute mass of the cells
        """
        self._mass = self.size_t * DataContainer().geometric.section * self.density.current_value

    @property
    def classical(self):
        """
        :return: boolean mask indicating which cells are classical
        """
        return self._classical

    @property
    def enriched(self):
        """
        :return: boolean mask indicating which cells are enriched
        """
        return ~self._classical

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

    @property
    def deviatoric_stress_new(self):
        """
        Return the new deviatoric stress tensor Sxx Syy Szz
        """
        return self._deviatoric_stress_new

    @property
    def equivalent_plastic_strain_rate(self):
        """
        Return the equivalent plastic strain rate
        """
        return self._equivalent_plastic_strain_rate

    @property
    def deviatoric_stress_current(self):
        """
        Return the current deviatoric stress tensor Sxx Syy Szz
        """
        return self._deviatoric_stress_current

    @property
    def plastic_strain_rate(self):
        """
        Return the current plastic stain rate tensor Dp_xx Dp_yy Dp_zz
        """
        return self._plastic_strain_rate

    def _compute_new_pressure_with_external_lib(self, density_current, density_new,
                                                pressure_current, pseudo_current,
                                                energy_current, energy_new, pressure_new, vson_new):
        """
        Computation of the set (internal energy, pressure, sound velocity) for v-e
        formulation thanks to external C library
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

    def compute_new_pressure(self, mask, dt):
        """
        Computation of the set (internal energy, pressure, sound velocity) for v-e formulation
        """
        elasticity_activated = np.logical_or(
            (DataContainer().material_target.constitutive_model.elasticity_model is not None),
            (DataContainer().material_projectile.constitutive_model.elasticity_model is not None))
        plasticity_activated = np.logical_or(
            (DataContainer().material_target.constitutive_model.plasticity_model is not None),
            (DataContainer().material_projectile.constitutive_model.plasticity_model is not None))

        if elasticity_activated or plasticity_activated:
            # si l'élasticité n'est pas activée, les grandeurs élastiques restent nulles.
            # L'opération est transparente pour les matériaux hydro.
            # Donc on peut ne pas distinguer les matériaux.
            self.energy.current_value[mask] += \
                OneDimensionCell.add_elastic_energy_method(dt,
                                                           self.density.current_value[mask],
                                                           self.density.new_value[mask],
                                                           self._deviatoric_stress_current[mask, :],
                                                           self._deviatoric_stress_new[mask, :],
                                                           self._deviatoric_strain_rate[mask, :])

        # Appel de l'équation d'état sur le projectile:
        mask_p = np.logical_and(mask, self.cell_in_projectile)
        if DataContainer().data_contains_a_projectile and mask_p.any():
            # Not sure there is cells in the intersection mask_classic / projectile
            # => test it before starting Newton
            self.energy.new_value[mask_p], self.pressure.new_value[mask_p], \
            self.sound_velocity.new_value[mask_p] = OneDimensionCell.apply_equation_of_state(
                self, DataContainer().material_projectile.constitutive_model.eos,
                self.density.current_value[mask_p], self.density.new_value[mask_p],
                self.pressure.current_value[mask_p], self.pressure.new_value[mask_p],
                self.energy.current_value[mask_p], self.energy.new_value[mask_p],
                self.pseudo.current_value[mask_p], self.sound_velocity.new_value[mask_p])

        # Appel de l'équation d'état sur la cible
        mask_t = np.logical_and(mask, self.cell_in_target)
        if DataContainer().data_contains_a_target and mask_t.any():
            # Not sure there is cells in the intersection mask_classic / target
            # => test it before starting Newton
            self.energy.new_value[mask_t], self.pressure.new_value[mask_t], \
            self.sound_velocity.new_value[mask_t] = OneDimensionCell.apply_equation_of_state(
                self, DataContainer().material_target.constitutive_model.eos,
                self.density.current_value[mask_t], self.density.new_value[mask_t],
                self.pressure.current_value[mask_t], self.pressure.new_value[mask_t],
                self.energy.current_value[mask_t], self.energy.new_value[mask_t],
                self.pseudo.current_value[mask_t],
                self.sound_velocity.new_value[mask_t])

    def compute_size(self, topologie, vecteur_coord_noeuds):
        """
        Computation of the cells initial length
        :param topologie : table to link nodes and cells index
        :param vecteur_coord_noeuds : array with nodes coordinates
        """
        connectivity = topologie.nodes_belonging_to_cell
        size = vecteur_coord_noeuds[connectivity[:, 1]] - vecteur_coord_noeuds[connectivity[:, 0]]
        cell_error = (size < 0)
        if cell_error.any():
            raise ValueError("La maille {:} a une longueur négative !".format(
                np.where(cell_error)[0]))
        self._size_t = size.flatten()

    def compute_new_size(self, topologie, vecteur_coord_noeuds, mask):
        """
        Computation of the cells length at time t+dt
        :param topologie : table to link nodes and cells index
        :param vecteur_coord_noeuds : array with the new nodes coordinates
        :param mask : array of boolean to identify classical cells
        """
        connectivity = topologie.nodes_belonging_to_cell
        size = vecteur_coord_noeuds[connectivity[:, 1]] - vecteur_coord_noeuds[connectivity[:, 0]]
        cell_error = (size < 0)
        if cell_error.any():
            raise ValueError("La maille {:} a une longueur négative !".format(
                np.where(cell_error)[0]))
        self._size_t_plus_dt[mask] = size[mask].flatten()

    def compute_new_density(self, mask):
        """
        Computation of the density of the cells at time t+dt using mass conservation principle
        :param mask : array of boolean to identify classical cells
        """
        self.density.new_value[mask] = self.density.current_value[mask] * \
                                       self.size_t[mask] / self.size_t_plus_dt[mask]

    def compute_new_pseudo(self, delta_t, mask):
        """
        Computation of cells artificial viscosity at time t+dt
        :param delta_t : time step
        :param mask : array of boolean to identify classical cells
        """
        self.pseudo.new_value[mask] = OneDimensionCell.compute_pseudo(
            delta_t, self.density.current_value[mask], self.density.new_value[mask],
            self.size_t_plus_dt[mask], self.sound_velocity.current_value[mask],
            DataContainer().numeric.a_pseudo, DataContainer().numeric.b_pseudo)

    def compute_new_time_step(self, mask):
        """
        Computation of the time step in the cells at time t+dt
        :param mask : array of boolean to identify classical cells
        """
        cfl = DataContainer().numeric.cfl
        cfl_pseudo = DataContainer().numeric.cfl_pseudo
        dt = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, self.density.current_value[mask],
                                                self.density.new_value[mask],
                                                self.size_t_plus_dt[mask],
                                                self.sound_velocity.new_value[mask],
                                                self.pseudo.current_value[mask],
                                                self.pseudo.new_value[mask])
        self._dt[mask] = dt

    def compute_shear_modulus(self):
        """
        Compute the shear modulus G according to the constitutive elasticity model in XDATA
        """
        pass

    def compute_yield_stress(self):
        """
        Compute the yield stress according to plasticity constitutve model in XDATA
        """
        pass

    def compute_complete_stress_tensor(self, mask):
        """
        Compute the Cauchy stress tensor (assemble pression et dï¿½viateur)
        :param mask : array of boolean to identify classical cells
        """
        for i in range(0, 3):
            self._stress[mask, i] = - (self.pressure.new_value[mask] + self.pseudo.new_value[mask])

        elasticity_activated = (
                DataContainer().material_target.constitutive_model.elasticity_model is not None or
                DataContainer().material_projectile.constitutive_model.elasticity_model is not None)
        plasticity_activated = (
                DataContainer().material_target.constitutive_model.plasticity_model is not None or
                DataContainer().material_projectile.constitutive_model.plasticity_model is not None)
        if elasticity_activated or plasticity_activated:
            self._stress[mask, :] += self._deviatoric_stress_new[mask, :]

    def compute_deviatoric_stress_tensor(self, mask, topologie, coord_noeud_new,
                                         vitesse_noeud_new, dt):
        """
        Compute the deviatoric part of the stress tensor
        :param mask : mask to select classical cells
        :param topologie : table of connectivity : link between cells and nodes id
        :param coord_noeud_new : array with new nodes coordinates
        :param vitesse_noeud_new : array with new nodes velocities
        :param dt : time step (staggered tn+1/2)
        """
        self.compute_shear_modulus()
        G = self.shear_modulus.current_value

        # Compute rotation rate tensor  and strain rate tensor: W = 0 en 1D
        self.compute_deviator_strain_rate(mask, dt, topologie, coord_noeud_new, vitesse_noeud_new)
        D = self._deviatoric_strain_rate

        # Rappel : S / dt * (-W * S + S * W) + 2. * G * deviateur_strain_rate[mask] * dt
        for i in range(0, 3):
            self._deviatoric_stress_new[mask, i] = \
                np.copy(self._deviatoric_stress_current[mask, i]) + 2. * G[mask] * D[mask, i] * dt

        # -----------------------------
        # pour être sur que la trace soit nulle
        trace = self._deviatoric_stress_new[mask, 0] + self._deviatoric_stress_new[mask, 1] \
                + self._deviatoric_stress_new[mask, 2]
        for i in range(0, 3):
            self._deviatoric_stress_new[mask, i] -= 1./3. * trace

    def compute_deviator_strain_rate(self, mask, dt, topologie, coord_noeud_new, vitesse_noeud_new):
        """
        Compute deviateur du taux de dï¿½formation
        :param mask : mask to select classical cells
        :param dt : time step
        :param topologie : table of connectivity : link between cells and nodes id
        :param coord_noeud_new : array with new nodes coordinates
        :param vitesse_noeud_new : array with new nodes velocities
        """
        connectivity = topologie.nodes_belonging_to_cell
        u_new = vitesse_noeud_new[connectivity][:, :, 0]
        # velocities of the left and right nodes à coté de cell
        x_new = coord_noeud_new[connectivity][:, :, 0]
        # coordinates of the left and right nodes à coté de cell

        # Calcul du déviateur de D
        self._deviatoric_strain_rate[mask, :] = \
            OneDimensionCell.general_method_deviator_strain_rate(mask, dt, x_new, u_new)

    def apply_plastic_corrector_on_deviatoric_stress_tensor(self, mask):
        """
        Correct the elastic trial of deviatoric stress tensor when plasticity criterion is activated
        :param mask : mask to identify the cells where plasticity should be applied
        (classical cells where plasticity criterion is activated)
        """
        invariant_J2_el = compute_J2(self.deviatoric_stress_new)
        # prédiction élastique avant le traitement de la plasticité
        radial_return = self.yield_stress.current_value / invariant_J2_el
        plasticity = radial_return < 1.
        plastic_mask = np.logical_and(mask, plasticity)
        for i in range(0, 3):
                self._deviatoric_stress_new[plastic_mask, i] *= radial_return[plastic_mask]

    def compute_plastic_strain_rate_tensor(self, mask, dt):
        """
        Compute the plastic strain rate tensor from elastic prediction and radial return
        (normal law for Von Mises plasticity)
        :param mask: mask to identify plastic cells
        :param dt: time step
        """
        # A faire avant apply_plastic_corrector_on_deviatoric_stress_tensor
        invariant_J2_el = compute_J2(self.deviatoric_stress_new)
        radial_return = self.yield_stress.current_value / invariant_J2_el
        plasticity = radial_return < 1.
        plastic_mask = np.logical_and(mask, plasticity)
        for i in range(0, 3):
            self._plastic_strain_rate[plastic_mask, i] = \
                (1 - radial_return[plastic_mask]) * self._deviatoric_stress_new[plastic_mask, i] / \
                (radial_return[plastic_mask] *
                 3 * self.shear_modulus.current_value[plastic_mask] * dt)


    def compute_equivalent_plastic_strain_rate(self, mask, dt):
        """
        Compute the plastic strain rate
        :param mask: array of bool to select cells of interest
        :param dt : float, time step staggered
        """
        invariant_J2_el = compute_J2(self.deviatoric_stress_new)
        # prédiction élastique avant le traitement de la plasticité
        G = self.shear_modulus.current_value
        plasticity = invariant_J2_el > self.yield_stress.current_value
        plastic_mask = np.logical_and(mask, plasticity)

        self._equivalent_plastic_strain_rate[plastic_mask] = \
            (invariant_J2_el[plastic_mask] - self.yield_stress.current_value[plastic_mask]) / \
            (3. * G[plastic_mask] * dt)

    def impose_pressure(self, ind_cell, pressure):
        """
        Pressure imposition
        :param ind_cell : index of the cells
        :param pressure : pressure value to be imposed
        """
        self.pressure.new_value[ind_cell] = pressure
        self._deviatoric_stress_new[ind_cell, :] = np.ones([3]) * pressure

    def increment_variables(self):
        """
        Increment cells variables from one iteration to another
        """
        super(OneDimensionCell,self).increment_variables()
        self._deviatoric_stress_current[:, :] = self._deviatoric_stress_new[:, :]
