# -*- coding: utf-8 -*-
"""
Implementation of the OneDimensionCell class
"""
import ctypes
import numpy as np

from xfv.src.cell import Cell
from xfv.src.solver.functionstosolve.vnrenergyevolutionforveformulation import (
    VnrEnergyEvolutionForVolumeEnergyFormulation)
from xfv.src.solver.newtonraphson import NewtonRaphson
from xfv.src.utilities.stress_invariants_calculation import compute_second_invariant

USE_INTERNAL_SOLVER = False
try:
    from launch_vnr_resolution_c import launch_vnr_resolution, MieGruneisenParams
except ImportError:
    USE_INTERNAL_SOLVER = True


# noinspection PyArgumentList
class OneDimensionCell(Cell):  # pylint: disable=too-many-public-methods
    """
    A class for one dimension cells
    """

    @classmethod
    def apply_equation_of_state(cls, cell, eos, density, density_new, pressure, pressure_new,
                                energy, energy_new, pseudo, cson_new):
        if not USE_INTERNAL_SOLVER:
            params = MieGruneisenParams(**eos.eos_param._asdict())
            pressure = pressure + 2. * pseudo
            launch_vnr_resolution(params, 1. / density, 1. / density_new, pressure, energy,
                                  energy_new, pressure_new, cson_new)
            return energy_new, pressure_new, cson_new
        else:
            my_variables = {'EquationOfState': eos,
                            'OldSpecificVolume': 1. / density,
                            'NewSpecificVolume': 1. / density_new,
                            'Pressure': (pressure + 2. * pseudo),
                            'OldEnergy': energy}
            cell._function_to_vanish.set_variables(my_variables)
            energy_new_value = cell._solver.compute_solution(energy)

            # Eos call to determine final pressure and sound speed values
            shape = energy_new.shape
            pressure_new_value = np.zeros(shape, dtype=np.float64, order='C')
            sound_velocity_new_value = np.zeros(shape, dtype=np.float64, order='C')
            dummy = np.zeros(shape, dtype=np.float64, order='C')
            my_variables['EquationOfState'].solve_volume_energy(
                my_variables['NewSpecificVolume'], energy_new_value, pressure_new_value,
                dummy, sound_velocity_new_value)

            if np.isnan(sound_velocity_new_value).any():
                negative_vson = np.where(np.isnan(sound_velocity_new_value))
                msg = "Sound speed square < 0 in cells {}\n".format(np.where(negative_vson))
                msg += "density = {}\n".format(my_variables['NewDensity'][negative_vson])
                msg += "energy = {}\n".format(energy_new_value[negative_vson])
                msg += "pressure = {}\n".format(pressure_new_value[negative_vson])
                raise ValueError(msg)

            cell._function_to_vanish.eraseVariables()
        return energy_new_value, pressure_new_value, sound_velocity_new_value

    @classmethod
    def add_elastic_energy_method(cls, dt,  # pylint: disable=invalid-name
                                  density_current, density_new,
                                  stress_dev_current, stress_dev_new, strain_rate_dev):
        """
        Take into account the additional term in internal energy due to elasticity
        """
        # The factor 1/2 from the mean between density_current et density_new simplifies with
        # the 1/2 coming from the mean between stress_dev_current et stress_dev_new
        energy_new_value = dt / (density_new + density_current)\
                           * ((stress_dev_new[:, 0]
                               + stress_dev_current[:, 0]) * strain_rate_dev[:, 0] +
                              (stress_dev_new[:, 1]
                               + stress_dev_current[:, 1]) * strain_rate_dev[:, 1] +
                              (stress_dev_new[:, 2]
                               + stress_dev_current[:, 2]) * strain_rate_dev[:, 2])
        return energy_new_value

    @classmethod
    def general_method_deviator_strain_rate(cls, mask, dt,  # pylint: disable=invalid-name
                                            x_new, u_new):
        """
        Compute the deviator of strain rate tensor (defined at the center of the cell)
        from the coordinates and velocities interpolated at the center of the cell
        :param mask : bool array to identify cells to be computed
        :param dt : time step (float)
        :param x_new : coordinates array at time n+1 of the left and right nodes of each cell
        Shape is array([coord_node_left, coord_node_right] * nbr_cells in the mask)
        :param u_new : velocity array at time n+1/2 of the left and right nodes of each cell
        Shape is array([velocity_node_left, velocity_node_right] * nbr_cells in the mask)
        x_new, u_new shape is (size(mask), 2)
        """
        strain_rate_dev = np.zeros([u_new.shape[0], 3])
        # Strain rate tensor
        x_demi = x_new - dt/2. * u_new
        D = (u_new[mask, 1] - u_new[mask, 0]) / (x_demi[mask, 1] - x_demi[mask, 0])  # Dxx
        # Cancel the trace to get the deviator part
        strain_rate_dev[mask, 0] = 2. / 3. * D
        strain_rate_dev[mask, 1] = - 1. / 3. * D
        strain_rate_dev[mask, 2] = - 1. / 3. * D
        return strain_rate_dev[mask]

    @classmethod
    def compute_pseudo(cls, delta_t: float, rho_old: np.array, rho_new: np.array,
                       size_new: np.array, cel_son: np.array, a_pseudo: float, b_pseudo: float):
        """
        Computation of artificial viscosity
        :param delta_t: time step
        :param rho_old: density at time t
        :param rho_new: density at time t+dt
        :param size_new: cell size at time t+dt
        :param cel_son: sound speed at time t
        :param a_pseudo: quadratic pseudo coefficient
        :param b_pseudo: linear pseudo coefficient
        """
        # pylint: disable=too-many-arguments
        # 8 arguments seems necessary
        vn = 1. / rho_old
        vnplusun = 1. / rho_new
        vnplusundemi = (vn + vnplusun) / 2.
        vpoint = (vnplusun - vn) / delta_t
        divu = vpoint / vnplusundemi
        pseudo = np.zeros(rho_old.shape, dtype=np.float64, order='C')
        mask = np.where(divu < 0.)
        pseudo[mask] = (1. / vnplusundemi[mask] *
                        (a_pseudo * size_new[mask] ** 2 * divu[mask] ** 2 +
                         b_pseudo * size_new[mask] * cel_son[mask] * abs(divu[mask])))
        return pseudo

    @classmethod
    def compute_time_step(cls, cfl, cfl_pseudo, rho_old, rho_new, size_new, sound_speed_new,
                          pseudo_old, pseudo_new):
        """
        Computation of the time step
        :param cfl : nombre cfl
        :param cfl_pseudo : cfl linked to the shock treatment stability condition
        :param rho_old : density at time t
        :param rho_new : density at time t+dt
        :param size_new : size of element
        :param sound_speed_new : sound velocity at time t+dt
        :param pseudo_old : artificial viscosity at time t
        :param pseudo_new : artificial viscosity at timet+dt
        """
        # pylint: disable=too-many-arguments
        local_cson = np.copy(sound_speed_new) ** 2
        mask_q = pseudo_new != 0.
        drho = np.abs((rho_new - rho_old) / rho_old)
        dpseudo = (pseudo_new - pseudo_old)
        mask_r = drho > 1.e-04
        mask_local_cson = np.logical_and(mask_q, mask_r)
        pseudo_sound_speed_square = np.abs(cfl_pseudo * dpseudo[mask_local_cson] /
                                           (rho_new[mask_local_cson] - rho_old[mask_local_cson]))
        local_cson[mask_local_cson] += pseudo_sound_speed_square
        local_cson **= 0.5
        delta_t = cfl * size_new / local_cson
        return delta_t

    def __init__(self, number_of_elements: int):
        """
        Build the array of cells (1D)
        :param number_of_elements: number of cells
        """
        super(OneDimensionCell, self).__init__(number_of_elements)

        # By default :all cells are classical (non enriched)
        self._classical = np.ones([number_of_elements, ], dtype=np.bool, order='C')
        self._enrichment_not_concerned = np.ones([number_of_elements, ], dtype=np.bool, order='C')

        self.dtc = []  # critical time step

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
        # EOS :
        self._target_eos = self.data.material_target.constitutive_model.eos.build_eos_obj()
        self._projectile_eos = None
        if self.data.data_contains_a_projectile:
            self._projectile_eos = \
                self.data.material_projectile.constitutive_model.eos.build_eos_obj()
        self._function_to_vanish = VnrEnergyEvolutionForVolumeEnergyFormulation()

        # Solver EOS
        self._solver = NewtonRaphson(self._function_to_vanish)

    def compute_mass(self):
        """
        Compute mass of the cells
        """
        self._mass = self.size_t * self.data.geometric.section * self.density.current_value

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

    def compute_new_pressure(self, mask, dt):  # pylint: disable=invalid-name
        """
        Computation of the set (internal energy, pressure, sound velocity) for v-e formulation
        """
        target_model = self.data.material_target.constitutive_model
        projectile_model = None
        if self.data.data_contains_a_projectile:
            projectile_model = self.data.material_projectile.constitutive_model
        elasticity_activated = np.logical_or(
            target_model.elasticity_model is not None,
            projectile_model and projectile_model.elasticity_model is not None)
        plasticity_activated = np.logical_or(
            target_model.plasticity_model is not None,
            projectile_model and projectile_model.plasticity_model is not None)

        if elasticity_activated or plasticity_activated:
            # if elasticity is not activated, then the elasticity fields remain null at each time.
            # Thus, no need to distinguish the projectile and target because operation is
            # transparent if no elasticity for one material
            self.energy.current_value[mask] += \
                OneDimensionCell.add_elastic_energy_method(dt,
                                                           self.density.current_value[mask],
                                                           self.density.new_value[mask],
                                                           self._deviatoric_stress_current[mask, :],
                                                           self._deviatoric_stress_new[mask, :],
                                                           self._deviatoric_strain_rate[mask, :])

        # Equation of state for the projectile:
        mask_p = np.logical_and(mask, self.cell_in_projectile)
        if self.data.data_contains_a_projectile and mask_p.any():
            # Not sure there is cells in the intersection mask_classic / projectile
            # => test it before starting Newton
            self.energy.new_value[mask_p], self.pressure.new_value[mask_p], \
            self.sound_velocity.new_value[mask_p] = OneDimensionCell.apply_equation_of_state(
                self, self._projectile_eos,
                self.density.current_value[mask_p], self.density.new_value[mask_p],
                self.pressure.current_value[mask_p], self.pressure.new_value[mask_p],
                self.energy.current_value[mask_p], self.energy.new_value[mask_p],
                self.pseudo.current_value[mask_p], self.sound_velocity.new_value[mask_p])

        # Equation of state for the target
        mask_t = np.logical_and(mask, self.cell_in_target)
        if self.data.data_contains_a_target and mask_t.any():
            # Not sure there is cells in the intersection mask_classic / target
            # => test it before starting Newton
            self.energy.new_value[mask_t], self.pressure.new_value[mask_t],\
            self.sound_velocity.new_value[mask_t] = OneDimensionCell.apply_equation_of_state(
                self, self._target_eos,
                self.density.current_value[mask_t], self.density.new_value[mask_t],
                self.pressure.current_value[mask_t], self.pressure.new_value[mask_t],
                self.energy.current_value[mask_t], self.energy.new_value[mask_t],
                self.pseudo.current_value[mask_t],
                self.sound_velocity.new_value[mask_t])

    def compute_size(self, topology, node_coord):
        """
        Computation of the cells initial length
        :param topology : table to link nodes and cells index
        :param node_coord : array with nodes coordinates
        """
        connectivity = topology.nodes_belonging_to_cell
        size = node_coord[connectivity[:, 1]] - node_coord[connectivity[:, 0]]
        cell_error = (size < 0)
        if cell_error.any():
            raise ValueError("Length of cell {:} is negative !".format(
                np.where(cell_error)[0]))
        self._size_t = size.flatten()

    def compute_new_size(self, topology, node_coord, mask):
        """
        Computation of the cells length at time t+dt
        :param topology : table to link nodes and cells index
        :param node_coord : array with the new nodes coordinates
        :param mask : array of boolean to identify classical cells
        """
        connectivity = topology.nodes_belonging_to_cell
        size = node_coord[connectivity[:, 1]] - node_coord[connectivity[:, 0]]
        cell_error = (size < 0)
        if cell_error[mask].any():
            raise ValueError("Length of cell {:} is negative !".format(
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
            self.data.numeric.a_pseudo, self.data.numeric.b_pseudo)

    def compute_new_time_step(self, mask):
        """
        Computation of the time step in the cells at time t+dt
        :param mask : array of boolean to identify classical cells
        """
        cfl = self.data.numeric.cfl
        cfl_pseudo = self.data.numeric.cfl_pseudo
        delta_t = OneDimensionCell.compute_time_step(cfl, cfl_pseudo,
                                                     self.density.current_value[mask],
                                                     self.density.new_value[mask],
                                                     self.size_t_plus_dt[mask],
                                                     self.sound_velocity.new_value[mask],
                                                     self.pseudo.current_value[mask],
                                                     self.pseudo.new_value[mask])
        self._dt[mask] = delta_t

    def compute_shear_modulus(self, shear_modulus_model, mask):
        """
        Compute the shear modulus G according to the constitutive elasticity model in XDATA
        :param shear_modulus_model : model to compute the shear modulus
        :param mask : mask to identify the cells to be computed
        """
        self.shear_modulus.new_value[mask] = \
            shear_modulus_model.compute(self.density.new_value[mask])

    def compute_yield_stress(self, yield_stress_model, mask):
        """
        Compute the yield stress according to plasticity constitutive model in XDATA
        :param yield_stress_model : model to compute the yield stress
        :param mask : mask to identify the cells to be computed
        """
        self.yield_stress.new_value[mask] = yield_stress_model.compute(self.density.new_value[mask])

    def compute_complete_stress_tensor(self):
        """
        Compute the Cauchy stress tensor (assemble pression et dï¿½viateur)
        """
        for i in range(0, 3):
            self._stress[:, i] = - (self.pressure.new_value + self.pseudo.new_value)
        self._stress += self._deviatoric_stress_new

    def compute_deviatoric_stress_tensor(self, mask, topology, node_coord_new,
                                         node_velocity_new, dt):  # pylint: disable=invalid-name
        """
        Compute the deviatoric part of the stress tensor
        :param mask : mask to select classical cells
        :param topology : table of connectivity : link between cells and nodes id
        :param node_coord_new : array with new nodes coordinates (time n+1)
        :param node_velocity_new : array with new nodes velocities (time n+1/2)
        :param dt : time step (staggered tn+1/2)
        """
        G = np.diag(self.shear_modulus.new_value[mask])  # pylint: disable=invalid-name

        # Compute rotation rate tensor and strain rate tensor: W = 0 en 1D
        self.compute_deviator_strain_rate(mask, dt, topology, node_coord_new, node_velocity_new)
        D = self._deviatoric_strain_rate  # pylint: disable=invalid-name

        # Rappel : S / dt * (-W * S + S * W) + 2. * G * deviateur_strain_rate[mask] * dt
        self._deviatoric_stress_new[mask] = self._deviatoric_stress_current[mask] + 2. * G.dot(D[mask]) * dt

        # -----------------------------
        # To ensure the trace to be null
        trace = np.sum(self._deviatoric_stress_new[mask], axis=1) / 3.
        full_trace = np.array([trace, trace, trace]).transpose()
        self._deviatoric_stress_new[mask] -= full_trace

    def compute_deviator_strain_rate(self, mask, dt,  # pylint: disable=invalid-name
                                     topology, node_coord_new, node_velocity_new):
        """
        Compute deviateur du taux de dï¿½formation
        :param mask : mask to select classical cells
        :param dt : time step
        :param topology : table of connectivity : link between cells and nodes id
        :param node_coord_new : array with new nodes coordinates
        :param node_velocity_new : array with new nodes velocities
        """
        connectivity = topology.nodes_belonging_to_cell
        u_new = node_velocity_new[connectivity][:, :, 0]
        # velocities of the left and right nodes belonging to cell
        x_new = node_coord_new[connectivity][:, :, 0]
        # coordinates of the left and right nodes belonging to cell

        # Compute the deviatoric strain rate tensor
        self._deviatoric_strain_rate[mask, :] = \
            OneDimensionCell.general_method_deviator_strain_rate(mask, dt, x_new, u_new)

    def apply_plastic_corrector_on_deviatoric_stress_tensor(self, mask):
        """
        Correct the elastic trial of deviatoric stress tensor when plasticity criterion is activated
        :param mask : mask to identify the cells where plasticity should be applied
        (classical cells where plasticity criterion is activated)
        """
        invariant_j2_el = compute_second_invariant(self.deviatoric_stress_new[mask, :])
        # elastic predictor before applying plasticity
        radial_return = self.yield_stress.new_value[mask] / invariant_j2_el
        for i in range(0, 3):
            self._deviatoric_stress_new[mask, i] *= radial_return

    def compute_plastic_strain_rate_tensor(self, mask, dt):  # pylint: disable=invalid-name
        """
        Compute the plastic strain rate tensor from elastic prediction and radial return
        (normal law for Von Mises plasticity)
        :param mask: mask to identify plastic cells
        :param dt: time step
        """
        # To be done before apply_plastic_corrector_on_deviatoric_stress_tensor
        invariant_j2_el = compute_second_invariant(self.deviatoric_stress_new[mask, :])
        radial_return = self.yield_stress.new_value[mask] / invariant_j2_el
        shear_modulus = self.shear_modulus.new_value[mask]
        for i in range(0, 3):
            self._plastic_strain_rate[mask, i] = \
                (1 - radial_return) * self._deviatoric_stress_new[mask, i] / \
                (radial_return * 3 * shear_modulus * dt)

    def compute_equivalent_plastic_strain_rate(self, mask, dt):  # pylint: disable=invalid-name
        """
        Compute the plastic strain rate
        :param mask: array of bool to select cells of interest
        :param dt : float, time step staggered
        """
        invariant_j2_el = compute_second_invariant(self.deviatoric_stress_new[mask, :])
        # elastic predictor before applying plasticity
        G = self.shear_modulus.new_value[mask]  # pylint: disable=invalid-name
        Y = self.yield_stress.new_value[mask]  # pylint: disable=invalid-name
        self._equivalent_plastic_strain_rate[mask] = (invariant_j2_el - Y) / (3. * G * dt)

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
        super(OneDimensionCell, self).increment_variables()
        self._deviatoric_stress_current[:, :] = self._deviatoric_stress_new[:, :]
