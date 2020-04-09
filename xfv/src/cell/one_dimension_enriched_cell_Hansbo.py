# -*- coding: utf-8 -*-
"""
Implementing the Element1dEnriched class for Hansbo&Hansbo enrichment
"""
import numpy as np

from xfv.src.cell.one_dimension_enriched_cell import OneDimensionEnrichedCell
from xfv.src.cell.one_dimension_cell import OneDimensionCell
from xfv.src.data.data_container import DataContainer
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.utilities.stress_invariants_calculation import compute_second_invariant


# noinspection PyArgumentList
class OneDimensionHansboEnrichedCell(OneDimensionEnrichedCell):
    """
    A collection of 1d enriched elements. Treatment for Hansbo enrichment
    """
    @classmethod
    def compute_discontinuity_borders_velocity(cls, disc, node_velocity):
        """
        Compute the velocitites of points at the discontinuity border
        :param disc: Discontinuity to be considered
        :param node_velocity : array with nodes velocity
        :return ug : velocity of the discontinuity left boundary
        :return ud : velocity of the discontinuity right boundary
        """
        epsilon = disc.position_in_ruptured_element
        u1g = node_velocity[disc.mask_in_nodes]
        u2d = node_velocity[disc.mask_out_nodes]
        u2g = disc.additional_dof_velocity_new[0]
        u1d = disc.additional_dof_velocity_new[1]

        correction_factor_gauche = 1.
        correction_factor_droite = 1.
        if DataContainer().material_target.failure_model.lump_mass_matrix == 'diag_eps':
            print("correction velocity discontinuity boundaries")
            correction_factor_gauche = 1. / epsilon
            correction_factor_droite = 1. / (1. - epsilon)
        u2g_correction = correction_factor_gauche * u2g
        u1d_correction = correction_factor_droite * u1d
        ug = u2g_correction * epsilon + u1g * (1. - epsilon)
        ud = u2d * epsilon + u1d_correction * (1. - epsilon)
        return ug, ud

    def __init__(self, number_of_elements):
        super(OneDimensionHansboEnrichedCell, self).__init__(number_of_elements)

    def initialize_additional_cell_dof(self, disc):
        """
        Values to intialize the right part fields when discontinuity disc is created
        """
        # Initialization of the current field value
        disc.additional_dof_density.current_value = \
            self.density.current_value[disc.mask_ruptured_cell]
        disc.additional_dof_pressure.current_value = \
            self.pressure.current_value[disc.mask_ruptured_cell]
        disc.additional_dof_sound_velocity.current_value = \
            self.sound_velocity.current_value[disc.mask_ruptured_cell]
        disc.additional_dof_energy.current_value = \
            self.energy.current_value[disc.mask_ruptured_cell]
        disc.additional_dof_artificial_viscosity.current_value = \
            self.pseudo.current_value[disc.mask_ruptured_cell]
        disc._additional_dof_deviatoric_stress_current = \
            self._deviatoric_stress_current[disc.mask_ruptured_cell]
        disc.additional_dof_shear_modulus.current_value = \
            self.shear_modulus.current_value[disc.mask_ruptured_cell]
        disc.additional_dof_yield_stress.current_value = \
            self.yield_stress.current_value[disc.mask_ruptured_cell]

        # Initialization of new value field
        # (so that the current value is not erased if the field is not updated in current step)
        disc.additional_dof_density.new_value = \
            self.density.new_value[disc.mask_ruptured_cell]
        disc.additional_dof_pressure.new_value = \
            self.pressure.new_value[disc.mask_ruptured_cell]
        disc.additional_dof_sound_velocity.new_value = \
            self.sound_velocity.new_value[disc.mask_ruptured_cell]
        disc.additional_dof_energy.new_value = \
            self.energy.new_value[disc.mask_ruptured_cell]
        disc.additional_dof_artificial_viscosity.new_value = \
            self.pseudo.new_value[disc.mask_ruptured_cell]
        disc.additional_dof_shear_modulus.new_value = \
            self.shear_modulus.new_value[disc.mask_ruptured_cell]
        disc.additional_dof_yield_stress.new_value = \
            self.yield_stress.new_value[disc.mask_ruptured_cell]
        disc._additional_dof_deviatoric_stress_new = \
            self._deviatoric_stress_new[disc.mask_ruptured_cell]
        # Other quantities initialization
        disc._additional_dof_deviatoric_strain_rate = \
            self._deviatoric_strain_rate[disc.mask_ruptured_cell]
        disc._additional_dof_stress = self._stress[disc.mask_ruptured_cell]
        disc._additional_dof_equivalent_plastic_strain_rate = \
            self._equivalent_plastic_strain_rate[disc.mask_ruptured_cell]

    @staticmethod
    def reconstruct_enriched_hydro_field(classical_field, enriched_field_name):
        """
        True field reconstruction from the classical and enriched fields
        :param classical_field: champ classique de type Field
        :param enriched_field_name: champ enrichi (str)
        :return: champ complet
        :rtype np.array
        """
        # To build the coordinates of cell field, the cracked cells of discontinuities must be
        # sorted by cell_id in order to manage shifts
        insertion_field = np.zeros([len(Discontinuity.discontinuity_list()), 2])
        # insertion_field is an array : ruptured_cell_id, right_field
        for disc in Discontinuity.discontinuity_list():
            enriched_field = getattr(disc, enriched_field_name)
            insertion_field[int(disc.label) - 1, 0] = int(disc.ruptured_cell_id)
            insertion_field[int(disc.label) - 1, 1] = enriched_field.current_value

        insertion_field = np.sort(insertion_field, 0)
        res = np.copy(classical_field.current_value)
        offset = 1
        for indice_cell_rompue in insertion_field[:, 0]:
            res = np.insert(res, int(indice_cell_rompue) + offset, insertion_field[offset - 1, 1])
            offset += 1
        return res

    @staticmethod
    def reconstruct_enriched_elasto_field(classical_field, enriched_field_name):
        """
        True field reconstruction from the classical and enriched fields
        :param classical_field: champ classique de type np.array
        :param enriched_field_name: champ enrichi (str)
        :return: champ complet
        :rtype np.array
        """
        # To build the coordinates of cell field, the cracked cells of discontinuities must be
        # sorted by cell_id in order to manage shifts
        insertion_field = np.zeros([len(Discontinuity.discontinuity_list()), 2])
        # insertion_field is an array ruptured_cell_id, right_field
        for disc in Discontinuity.discontinuity_list():
            enriched_field = getattr(disc, enriched_field_name)
            insertion_field[int(disc.label) - 1, 0] = int(disc.ruptured_cell_id)
            insertion_field[int(disc.label) - 1, 1] = enriched_field[0][0]
        insertion_field = np.sort(insertion_field, 0)
        res = np.copy(classical_field[:, 0])
        offset = 1
        for indice_cell_rompue in insertion_field[:, 0]:
            res = np.insert(res, int(indice_cell_rompue) + offset, insertion_field[offset - 1, 1])
            offset += 1
        return res

    @property
    def pressure_field(self):
        """
        :return: pressure field
        :rtype: np.array
        """
        return self.reconstruct_enriched_hydro_field(self.pressure, "additional_dof_pressure")

    @property
    def density_field(self):
        """
        :return: density field
        :rtype: np.array
        """
        return self.reconstruct_enriched_hydro_field(self.density, "additional_dof_density")

    @property
    def energy_field(self):
        """
        :return: energy field
        :rtype: np.array
        """
        return self.reconstruct_enriched_hydro_field(self.energy, "additional_dof_energy")

    @property
    def artificial_viscosity_field(self):
        """
        :return: artificial viscosity field
        :rtype: np.array
        """
        return self.reconstruct_enriched_hydro_field(self.pseudo,
                                                     "additional_dof_artificial_viscosity")

    @property
    def stress_xx_field(self):
        """
        :return: sigma_xx field
        :rtype: np.array
        """
        return self.reconstruct_enriched_elasto_field(self.stress, "additional_dof_stress")

    @property
    def deviatoric_stress_field(self):
        """
        :return: (deviateur de sigma)_xx field
        :rtype: np.array
        """
        return self.reconstruct_enriched_elasto_field(self.deviatoric_stress_current,
                                                      "additional_dof_deviatoric_stress_current")

    def compute_enriched_elements_new_pressure(self, delta_t):
        """
        Compute pressure, internal energy and sound velocity in left and right parts of
        the enriched elements
        :param delta_t : time step
        """
        target_model = self.data.material_target.constitutive_model
        projectile_model = None
        if self.data.data_contains_a_projectile:
            projectile_model = self.data.material_projectile.constitutive_model
        elasticity_activated = \
            np.logical_or(target_model.elasticity_model is not None,
                          projectile_model and projectile_model.elasticity_model is not None)
        plasticity_activated = \
            np.logical_or(target_model.plasticity_model is not None,
                          projectile_model and projectile_model.plasticity_model is not None)

        for disc in Discontinuity.discontinuity_list():
            mask = disc.mask_ruptured_cell
            if elasticity_activated or plasticity_activated:
                self.energy.current_value[mask] += \
                    OneDimensionCell.add_elastic_energy_method(delta_t,
                                                               self.density.current_value[mask],
                                                               self.density.new_value[mask],
                                                               self.deviatoric_stress_current[mask],
                                                               self.deviatoric_stress_new[mask],
                                                               self._deviatoric_strain_rate[mask])
                disc.additional_dof_energy.current_value += \
                    OneDimensionCell.add_elastic_energy_method(
                        delta_t, disc.additional_dof_density.current_value,
                        disc.additional_dof_density.new_value,
                        disc.additional_dof_deviatoric_stress_current,
                        disc.additional_dof_deviatoric_stress_new,
                        disc.additional_dof_deviatoric_strain_rate)

            # Initialize local parameters :
            density_left = self.density.current_value[mask]
            density_left_new = self.density.new_value[mask]
            pressure_left = self.pressure.current_value[mask]
            pressure_left_new = self.pressure.new_value[mask]
            energy_left = self.energy.current_value[mask]
            energy_left_new = self.energy.new_value[mask]
            pseudo_left = self.pseudo.current_value[mask]
            cson_left_new = self.sound_velocity.new_value[mask]
            #
            density_right = disc.additional_dof_density.current_value
            density_right_new = disc.additional_dof_density.new_value
            pressure_right = disc.additional_dof_pressure.current_value
            pressure_right_new = disc.additional_dof_pressure.new_value
            energy_right = disc.additional_dof_energy.current_value
            energy_right_new = disc.additional_dof_energy.new_value
            pseudo_right = disc.additional_dof_artificial_viscosity.current_value
            cson_right_new = disc.additional_dof_sound_velocity.new_value
            # Call EOS :
            energy_new_left_value, pressure_new_left_value, sound_velocity_new_left_value = \
                OneDimensionCell.apply_equation_of_state(
                    self, self._target_eos,
                    density_left, density_left_new, pressure_left,
                    pressure_left_new, energy_left, energy_left_new,
                    pseudo_left, cson_left_new)

            energy_new_right_value, pressure_new_right_value, sound_velocity_new_right_value = \
                OneDimensionCell.apply_equation_of_state(
                    self, self._target_eos,
                    density_right, density_right_new, pressure_right,
                    pressure_right_new, energy_right, energy_right_new,
                    pseudo_right, cson_right_new)

            # Save results :
            self.pressure.new_value[mask] = pressure_new_left_value
            disc.additional_dof_pressure.new_value = pressure_new_right_value
            self.energy.new_value[mask] = energy_new_left_value
            disc.additional_dof_energy.new_value = energy_new_right_value
            self.sound_velocity.new_value[mask] = sound_velocity_new_left_value
            disc.additional_dof_sound_velocity.new_value = sound_velocity_new_right_value

    def compute_enriched_elements_new_part_size(self, time_step, node_velocity):
        """
        Computethe new size of each ruptured element part (left size and right size)
        :param time_step: time step
        :param node_velocity: array, node velocities
        """
        for disc in Discontinuity.discontinuity_list():
            ug, ud = OneDimensionHansboEnrichedCell.compute_discontinuity_borders_velocity(
                disc, node_velocity)
            u1g = node_velocity[disc.mask_in_nodes]
            u2d = node_velocity[disc.mask_out_nodes]
            OneDimensionEnrichedCell.compute_new_left_right_size(time_step, disc, u1g, u2d, ug, ud)

    def compute_enriched_elements_new_density(self):
        """
        Compute the new densities for left and right parts of the ruptured element
        (from mass conservation equation)
        """
        for disc in Discontinuity.discontinuity_list():
            mask = disc.ruptured_cell_id
            density_left = self.density.current_value[mask]
            density_right = disc.additional_dof_density.current_value
            self.density.new_value[mask], disc.additional_dof_density.new_value = \
                OneDimensionEnrichedCell.compute_new_left_right_density(density_left,
                                                                        density_right, disc)

    def compute_enriched_elements_new_pseudo(self, delta_t):
        """
        Calcule les nouvelles pseudo viscositï¿½s gauche et droite pour les ï¿½lï¿½ments enrichis
        ï¿½ partir de la methode compute_new_pseudo de OneDimensionCell avec les nouvelles
        valeurs enrichies
        :param delta_t: time_step
        """
        for disc in Discontinuity.discontinuity_list():
            mask_in = disc.ruptured_cell_id
            # Partie gauche :
            density_left = np.array([self.density.current_value[mask_in]])
            density_left_new = np.array([self.density.new_value[mask_in]])
            sound_velocity_left = np.array([self.sound_velocity.current_value[mask_in]])
            pseudo_left_new = OneDimensionCell.compute_pseudo(
                delta_t, density_left, density_left_new, disc.left_part_size.new_value,
                sound_velocity_left, self.data.numeric.a_pseudo,
                self.data.numeric.b_pseudo)
            # Partie droite :
            density_right = disc.additional_dof_density.current_value
            density_right_new = disc.additional_dof_density.new_value
            sound_velocity_right = disc.additional_dof_sound_velocity.current_value
            pseudo_right_new = OneDimensionCell.compute_pseudo(
                delta_t, density_right, density_right_new,
                disc.right_part_size.new_value, sound_velocity_right,
                self.data.numeric.a_pseudo, self.data.numeric.b_pseudo)
            self.pseudo.new_value[mask_in] = pseudo_left_new
            disc.additional_dof_artificial_viscosity.new_value = pseudo_right_new

    def compute_enriched_stress_tensor(self):
        """
        Compute the complete enriched stress tensor : sigma = -(p+q) I + S
        """
        for disc in Discontinuity.discontinuity_list():
            mask = disc.ruptured_cell_id
            self._stress[mask, 0] = self._deviatoric_stress_new[mask, 0] - \
                                    (self.pressure.new_value[mask] + self.pseudo.new_value[mask])
            self._stress[mask, 1] = self._deviatoric_stress_new[mask, 1] - \
                                    (self.pressure.new_value[mask] + self.pseudo.new_value[mask])
            self._stress[mask, 2] = self._deviatoric_stress_new[mask, 2] - \
                                    (self.pressure.new_value[mask] + self.pseudo.new_value[mask])

            # index 0 after pressure and pseudo viscosity fields to get scalar values
            # Each line is a scalar equation
            disc.additional_dof_stress[:, 0] = \
                disc.additional_dof_deviatoric_stress_new[:, 0] - \
                (disc.additional_dof_pressure.new_value[0] +
                 disc.additional_dof_artificial_viscosity.new_value[0])
            disc.additional_dof_stress[:, 1] = \
                disc.additional_dof_deviatoric_stress_new[:, 1] - \
                (disc.additional_dof_pressure.new_value[0] +
                 disc.additional_dof_artificial_viscosity.new_value[0])
            disc.additional_dof_stress[:, 2] = \
                disc.additional_dof_deviatoric_stress_new[:, 2] - \
                (disc.additional_dof_pressure.new_value[0] +
                 disc.additional_dof_artificial_viscosity.new_value[0])

    def compute_enriched_deviatoric_strain_rate(self, dt, node_coord_new, node_velocity_new):
        """
        Compute devaiateur du taux de dï¿½formation
        :param dt : time step
        :param node_coord_new : array, new nodes coordinates
        :param node_velocity_new : array, new nodes velocity
        """
        for disc in Discontinuity.discontinuity_list():
            mask_nodes = disc.mask_in_nodes + disc.mask_out_nodes
            mask_cells = disc.ruptured_cell_id

            u_discg_new, u_discd_new = \
                OneDimensionHansboEnrichedCell.compute_discontinuity_borders_velocity(
                    disc, node_velocity_new)
            u_noeuds_new = node_velocity_new[mask_nodes]  # left / right node velocity at time n+1
            x_noeuds_new = node_coord_new[mask_nodes]  # left / right node coordinates at time n+1

            # Creation of structure left - right data to call general_method_deviator_strain_rate
            # Left part cell : node_g - left boundary of discontinuity
            xg_new = np.array([x_noeuds_new[0],
                               x_noeuds_new[0] + disc.left_part_size.new_value]).reshape(1, 2)
            ug_new = np.array([u_noeuds_new[0], u_discg_new]).reshape(1, 2)
            # Compute the deviatoric strain rate tensor for left part
            self._deviatoric_strain_rate[mask_cells] = \
                OneDimensionCell.general_method_deviator_strain_rate(
                np.array([True]), dt, xg_new, ug_new)  # np.array(True) to be consistent

            # Creation of structure left - right data to call general_method_deviator_strain_rate
            # Left part cell : right boundary of discontinuity - node_right
            xd_new = np.array([x_noeuds_new[1] - disc.right_part_size.new_value,
                               x_noeuds_new[1]]).reshape(1, 2)
            ud_new = np.array([u_discd_new, u_noeuds_new[1]]).reshape(1, 2)
            # Compute the deviatoric strain rate tensor for right part
            disc._additional_dof_deviatoric_strain_rate = \
                OneDimensionCell.general_method_deviator_strain_rate(
                np.array([True]), dt, xd_new, ud_new)  # np.array(True) to be consistent

    def compute_enriched_deviatoric_stress_tensor(self, node_coord_new, node_velocity_new, dt):
        """
        Compute the deviatoric part of the stress tensor
        :param node_coord_new : array, new nodes coordinates
        :param node_velocity_new : array, new nodes velocity
        :param dt : float, time step
        """

        self.compute_enriched_deviatoric_strain_rate(dt, node_coord_new, node_velocity_new)

        # Compute rotation rate tensor : W = 0 en 1D

        self.compute_enriched_shear_modulus()

        for disc in Discontinuity.discontinuity_list():
            # Left part
            mask_in = disc.ruptured_cell_id
            G = self.shear_modulus.new_value[mask_in]
            self._deviatoric_stress_new[mask_in, 0] = \
                self._deviatoric_stress_current[mask_in, 0] + \
                2. * G * self._deviatoric_strain_rate[mask_in, 0] * dt
            self._deviatoric_stress_new[mask_in, 1] = \
                self._deviatoric_stress_current[mask_in, 1] + \
                2. * G * self._deviatoric_strain_rate[mask_in, 1] * dt
            self._deviatoric_stress_new[mask_in, 2] = \
                self._deviatoric_stress_current[mask_in, 2] + \
                2. * G * self._deviatoric_strain_rate[mask_in, 2] * dt
            # Right part
            Gd = disc.additional_dof_shear_modulus.new_value
            disc._additional_dof_deviatoric_stress_new[:, 0] = \
                disc.additional_dof_deviatoric_stress_current[:, 0] + \
                2. * Gd * disc.additional_dof_deviatoric_strain_rate[:, 0] * dt
            disc._additional_dof_deviatoric_stress_new[:, 1] = \
                disc.additional_dof_deviatoric_stress_current[:, 1] + \
                2. * Gd * disc.additional_dof_deviatoric_strain_rate[:, 1] * dt
            disc._additional_dof_deviatoric_stress_new[:, 2] = \
                disc.additional_dof_deviatoric_stress_current[:, 2] +\
                2. * Gd * disc.additional_dof_deviatoric_strain_rate[:, 2] * dt

    def compute_enriched_shear_modulus(self):
        """
        Compute the shear modulus for ruptured cell
        """
        # TODO : interroger le package rheology
        for disc in Discontinuity.discontinuity_list():
            mask = disc.ruptured_cell_id
            self.shear_modulus.new_value[mask] = \
                self.data.material_target.initial_values.shear_modulus_init
            disc.additional_dof_shear_modulus.new_value = \
                self.data.material_target.initial_values.shear_modulus_init

    def apply_plastic_correction_on_enriched_deviatoric_stress_tensor(self, mask):
        """
        Correct the elastic trial of deviatoric stress tensor when plasticity criterion is activated
        :param mask: mask to select plastic (enriched) cells where plasticity should be applied
        """
        for disc in Discontinuity.discontinuity_list():
            # Left part of the cracked cell :
            mask = np.logical_and(disc.mask_ruptured_cell, mask)
            invariant_j2_el = compute_second_invariant(self.deviatoric_stress_new)
            radial_return = self.yield_stress.current_value[mask] / invariant_j2_el[mask]
            self._deviatoric_stress_new[mask, 0] *= radial_return
            self._deviatoric_stress_new[mask, 1] *= radial_return
            self._deviatoric_stress_new[mask, 2] *= radial_return

            # Right part of the cracked cell :
            if disc.plastic_cells:  # if right part is plastic
                invariant_j2_el_right = \
                    compute_second_invariant(disc.additional_dof_deviatoric_stress_new)
                radial_return = \
                    disc.additional_dof_yield_stress.current_value / invariant_j2_el_right
                disc._additional_dof_deviatoric_stress_new[:, 0] *= radial_return
                disc._additional_dof_deviatoric_stress_new[:, 1] *= radial_return
                disc._additional_dof_deviatoric_stress_new[:, 2] *= radial_return

    def compute_enriched_equivalent_plastic_strain_rate(self, mask, dt):
        """
        Compute the plastic strain rate
        :param dt : time step
        :param mask: mask to select plastic cells where plasticity should be applied
        """
        for disc in Discontinuity.discontinuity_list():
            # Left part
            mask = np.logical_and(disc.mask_ruptured_cell, mask)
            invariant_j2_el = \
                compute_second_invariant(self._deviatoric_stress_new)  # elastic predictor
            G = self.shear_modulus.current_value[mask]  # pylint disable=invalid-name
            self._equivalent_plastic_strain_rate[mask] += \
                (invariant_j2_el[mask] - self.yield_stress.current_value[mask]) / (3. * G * dt)
            # Right part :
            if disc.plastic_cells:  # if the cracked cell is plastic :
                invariant_J2_el_right = \
                    compute_second_invariant(
                        disc.additional_dof_deviatoric_stress_new)  # elastic predictor
                Gd = disc.additional_dof_shear_modulus.current_value  # pylint disable=invalid-name
                disc._additional_dof_equivalent_plastic_strain_rate += \
                    (invariant_J2_el_right -
                     disc.additional_dof_yield_stress.current_value) / (3. * Gd * dt)

    def compute_enriched_plastic_strain_rate(self, mask, dt):
        """
        Compute the plastic strain rate tensor from elastic prediction and radial return
        (normal law for Von Mises plasticity) in cracked cells
        :param mask: mask to identify plastic cells
        :param dt: time step
        """
        for disc in Discontinuity.discontinuity_list():
            mask = np.logical_and(disc.mask_ruptured_cell, mask)
            invariant_j2_el = compute_second_invariant(self.deviatoric_stress_new)
            radial_return = self.yield_stress.current_value[mask] / invariant_j2_el[mask]
            for i in range(0, 3):
                self._plastic_strain_rate[mask, i] = \
                    (1 - radial_return) / (radial_return *
                                           3 * self.shear_modulus.current_value[mask] *
                                           dt) * self._deviatoric_stress_new[mask, i]

            # Partie droite :
            if disc.plastic_cells:  # si la cell rompue qu'on regarde est plastique :
                invariant_j2_el_right = \
                    compute_second_invariant(disc.additional_dof_deviatoric_stress_new)
                radial_return = \
                    disc.additional_dof_yield_stress.current_value / invariant_j2_el_right
                for i in range(0, 3):
                    disc._additional_dof_plastic_strain_rate[:, i] = \
                        (1- radial_return) / (radial_return * 3 *
                                              disc.additional_dof_shear_modulus.current_value *
                                              dt) * disc._additional_dof_deviatoric_stress_new[:, i]

    def compute_enriched_yield_stress(self):
        """
        Compute the yield stress for ruptured cells
        """
        for disc in Discontinuity.discontinuity_list():
            mask = disc.ruptured_cell_id
            # TODO : interroger le package rheology
            self.yield_stress.new_value[mask] = \
                self.data.material_target.initial_values.yield_stress_init
            disc.additional_dof_yield_stress.new_value = \
                self.data.material_target.initial_values.yield_stress_init

    def compute_enriched_elements_new_time_step(self):
        """
        Compute the new local time step.
        The calculation is equivalent to a remeshing time step and thus underestimates the
        time step for the enriched cells
        """
        cfl = self.data.numeric.cfl
        cfl_pseudo = self.data.numeric.cfl_pseudo

        for disc in Discontinuity.discontinuity_list():
            mask_in = disc.ruptured_cell_id

            # Left part
            # To be compatible with the method in Cell1D : take arrays of size 1
            density_left = np.array([self.density.current_value[mask_in]])
            density_left_new = np.array([self.density.new_value[mask_in]])
            sound_velocity_left_new = np.array([self.sound_velocity.new_value[mask_in]])
            pseudo_left = np.array([self.pseudo.current_value[mask_in]])
            pseudo_left_new = np.array([self.pseudo.new_value[mask_in]])

            dt_g = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, density_left,
                                                      density_left_new,
                                                      disc.left_part_size.new_value,
                                                      sound_velocity_left_new,
                                                      pseudo_left, pseudo_left_new)

            # Right part
            density_right = disc.additional_dof_density.current_value
            density_right_new = disc.additional_dof_density.new_value
            sound_velocity_right_new = disc.additional_dof_sound_velocity.new_value
            pseudo_right = disc.additional_dof_artificial_viscosity.current_value
            pseudo_right_new = disc.additional_dof_artificial_viscosity.new_value

            dt_d = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, density_right,
                                                      density_right_new,
                                                      disc.right_part_size.new_value,
                                                      sound_velocity_right_new,
                                                      pseudo_right, pseudo_right_new)

            self._dt[mask_in] = min(dt_g, dt_d)
