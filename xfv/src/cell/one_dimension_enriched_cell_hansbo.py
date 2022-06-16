# -*- coding: utf-8 -*-
"""
Implementing the Element1dEnriched class for Hansbo&Hansbo enrichment
"""
import os
from typing import Tuple

import numpy as np

from xfv.src.cell.one_dimension_cell import OneDimensionCell, Cell
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.utilities.stress_invariants_calculation import compute_second_invariant
from xfv.src.fields.field import Field


# noinspection PyArgumentList
class OneDimensionHansboEnrichedCell(OneDimensionCell):  # pylint: disable=too-many-public-methods
    """
    A collection of 1d enriched elements. Treatment for Hansbo enrichment
    """
    @staticmethod
    def _compute_discontinuity_borders_velocity(epsilon: np.array,
                                                u1g: np.array,
                                                u1d: np.array,
                                                u2g: np.array,
                                                u2d: np.array) -> Tuple[np.array, np.array]:
        """
        Compute the velocities of points at the discontinuity border

        :param epsilon: relative position of the discontinuity inside the cell
        :param u1g: classic velocity on left node (inside node)
        :param u1d: additional dof velocity on left node
        :param u2g: additional dof velocity on right node
        :param u2d: classical velocity on right node (outside)
        :return ug: velocity of the discontinuity left boundary
        :return ud: velocity of the discontinuity right boundary
        """
        ug = u2g * epsilon + u1g * (1. - epsilon)  # pylint: disable=invalid-name
        ud = u2d * epsilon + u1d * (1. - epsilon)  # pylint: disable=invalid-name
        return ug, ud

    @classmethod
    def compute_discontinuity_borders_velocity(cls, disc, node_velocity):
        """
        Compute the velocities of points at the discontinuity border

        :param disc: Discontinuity to be considered
        :param node_velocity: array with nodes velocity
        :return ug: velocity of the discontinuity left boundary
        :return ud: velocity of the discontinuity right boundary
        """
        epsilon = disc.position_in_ruptured_element
        u1g = node_velocity[disc.in_nodes]
        u2d = node_velocity[disc.out_nodes]
        u2g = disc.enr_velocity_new[0]
        u1d = disc.enr_velocity_new[1]
        # ug, ud = cls._compute_discontinuity_borders_velocity(epsilon, u1g, u1d, u2g, u2d)
        ug = u2g * epsilon + u1g * (1. - epsilon)  # pylint: disable=invalid-name
        ud = u2d * epsilon + u1d * (1. - epsilon)  # pylint: disable=invalid-name
        return ug, ud

    def __init__(self, n_cells: int):
        """
        Build the class OneDimensionHansboEnrichedCell

        :param n_cells: total number of cells
        """
        super().__init__(n_cells)
        #
        print(self._fields_manager)
        self._classical = np.ones(n_cells, dtype=bool, order='C')

        # Cell parts geometry
        self._left_part_size = Field(n_cells, np.zeros([n_cells]), np.zeros([n_cells]))
        self._right_part_size = Field(n_cells, np.zeros([n_cells]), np.zeros([n_cells]))
        self._enr_coordinates_x = np.zeros([n_cells], dtype=float)

        # Additional dof: thermodynamics
        self._enr_pressure = Field(n_cells, np.zeros([n_cells]), np.zeros([n_cells]))
        self._enr_density = Field(n_cells, np.zeros([n_cells]), np.zeros([n_cells]))
        self._enr_energy = Field(n_cells, np.zeros([n_cells]), np.zeros([n_cells]))
        self._enr_artificial_viscosity = Field(n_cells, np.zeros([n_cells]), np.zeros([n_cells]))
        self._enr_sound_velocity = Field(n_cells, np.zeros([n_cells]), np.zeros([n_cells]))

        # Additional dof : elasticity / plasticity
        self._enr_shear_modulus = Field(n_cells, np.zeros([n_cells]), np.zeros([n_cells]))
        self._enr_yield_stress = Field(n_cells, np.zeros([n_cells]), np.zeros([n_cells]))
        self._enr_stress = np.zeros([n_cells, 3])
        self._enr_deviatoric_stress_current = np.zeros([n_cells, 3])
        self._enr_deviatoric_stress_new = np.zeros([n_cells, 3])
        self._enr_deviatoric_strain_rate = np.zeros([n_cells, 3])
        self._enr_equivalent_plastic_strain_rate = np.zeros([n_cells])
        self._enr_plastic_strain_rate = np.zeros([n_cells, 3])

        # Indicator right part of enr cell is plastic
        self.plastic_enr_cells = np.zeros([n_cells], dtype=bool)

    def initialize_additional_cell_dof(self, disc: Discontinuity):
        """
        Values to initialize the right part fields when discontinuity disc is created

        :param disc: the current discontinuity
        """
        enr_cell = disc.get_ruptured_cell_id
        # Initialization of the current field value
        self.enr_density.current_value[enr_cell] = self.density.current_value[enr_cell]
        self.enr_pressure.current_value[enr_cell] = self.pressure.current_value[enr_cell]
        self.enr_sound_velocity.current_value[enr_cell] = self.sound_velocity.current_value[enr_cell]
        self.enr_energy.current_value[enr_cell] = self.energy.current_value[enr_cell]
        self.enr_artificial_viscosity.current_value[enr_cell] = self.pseudo.current_value[enr_cell]
        self._enr_deviatoric_stress_current[enr_cell] = self._deviatoric_stress_current[enr_cell]
        self.enr_shear_modulus.current_value[enr_cell] = self.shear_modulus.current_value[enr_cell]
        self.enr_yield_stress.current_value[enr_cell] = self.yield_stress.current_value[enr_cell]
        self._enr_coordinates_x[enr_cell] = self._coordinates_x[enr_cell]

        # Initialization of new value field
        # (so that the current value is not erased if the field is not updated in current step)
        self.enr_density.new_value[enr_cell] = self.density.new_value[enr_cell]
        self.enr_pressure.new_value[enr_cell] = self.pressure.new_value[enr_cell]
        self.enr_sound_velocity.new_value[enr_cell] = self.sound_velocity.new_value[enr_cell]
        self.enr_energy.new_value[enr_cell] = self.energy.new_value[enr_cell]
        self.enr_artificial_viscosity.new_value[enr_cell] = self.pseudo.new_value[enr_cell]
        self.enr_shear_modulus.new_value[enr_cell] = self.shear_modulus.new_value[enr_cell]
        self.enr_yield_stress.new_value[enr_cell] = self.yield_stress.new_value[enr_cell]
        self._enr_deviatoric_stress_new[enr_cell] = self._deviatoric_stress_new[enr_cell]
        # Other quantities initialization
        self._enr_deviatoric_strain_rate[enr_cell] = self._deviatoric_strain_rate[enr_cell]
        self._enr_stress[enr_cell] = self._stress[enr_cell]
        self._enr_equivalent_plastic_strain_rate[enr_cell] = self._equivalent_plastic_strain_rate[enr_cell]

    def reconstruct_enriched_hydro_field(self, classical_field: Field, enriched_field_name: str):
        """
        True field reconstruction from the classical and enriched fields

        :param classical_field: classical field
        :param enriched_field_name: name of the enriched field
        :return: complete field
        :rtype: np.array
        """
        # To build the coordinates of cell field, the cracked cells of discontinuities must be
        # sorted by cell_id in order to manage shifts
        enriched_field = getattr(self, enriched_field_name)
        insertion_field = np.zeros([len(Discontinuity.discontinuity_list()), 2])
        insertion_field[:, 0] = np.arange(self.number_of_cells)[self.enriched]  # filtered indexes
        insertion_field[:, 1] = enriched_field.current_value[self.enriched]  # filtered field
        res = np.copy(classical_field.current_value)
        offset = 1
        for indice_cell_rompue in insertion_field[:, 0]:
            res = np.insert(res, int(indice_cell_rompue) + offset, insertion_field[offset - 1, 1])
            offset += 1
        return res

    def reconstruct_enriched_elasto_field(self, classical_field: np.array,
                                          enriched_field_name: str):
        """
        True field reconstruction from the classical and enriched fields

        :param classical_field: classical field
        :param enriched_field_name: name of the enriched field
        :return: complete field
        :rtype: np.array
        """
        # To build the coordinates of cell field, the cracked cells of discontinuities must be
        # sorted by cell_id in order to manage shifts
        enriched_field = getattr(self, enriched_field_name)
        insertion_field = np.zeros([len(Discontinuity.discontinuity_list()), 2])
        insertion_field[:, 0] = np.arange(self.number_of_cells)[self.enriched]  # filtered indexes
        insertion_field[:, 1] = enriched_field[self.enriched, 0]  # filtered field_xx
        res = np.copy(classical_field[:, 0])
        offset = 1
        for indice_cell_rompue in insertion_field[:, 0]:
            res = np.insert(res, int(indice_cell_rompue) + offset, insertion_field[offset - 1, 1])
            offset += 1
        return res

    @property
    def left_part_size(self):
        """
        Accessor on the size of the left part of cracked cell field
        """
        return self._left_part_size

    @property
    def right_part_size(self):
        """
        Accessor on the size of the right part of cracked cell field
        """
        return self._right_part_size

    @property
    def enr_pressure(self):
        """
        Accessor on the right part of cracked cell pressure field
        """
        return self._enr_pressure

    @property
    def enr_density(self):
        """
        Accessor on the right part of cracked cell density field
        """
        return self._enr_density

    @property
    def enr_sound_velocity(self):
        """
        Accessor on the right part of cracked cell sound speed field
        """
        return self._enr_sound_velocity

    @property
    def enr_energy(self):
        """
        Accessor on the right part of cracked cell internal energy field
        """
        return self._enr_energy

    @property
    def enr_artificial_viscosity(self):
        """
        Accessor on the right part of cracked cell artificial viscosity field
        """
        return self._enr_artificial_viscosity

    @property
    def enr_stress(self):
        """
        Accessor on the right part of cracked cell stress at time t
        """
        return self._enr_stress

    @property
    def enr_stress_xx(self):
        """
        Accessor on the right part of cracked cell stress at time t
        """
        return self._enr_stress[:, 0]

    @property
    def enr_deviatoric_stress_current(self):
        """
        Accessor on the right part of cracked cell deviatoric stress at time t
        """
        return self._enr_deviatoric_stress_current

    @property
    def enr_deviatoric_stress_new(self):
        """
        Accessor on the right part of cracked cell deviatoric stress at time t+dt
        """
        return self._enr_deviatoric_stress_new

    @property
    def enr_deviatoric_strain_rate(self):
        """
        Accessor on the right part of cracked cell deviatoric strain rate at time t
        """
        return self._enr_deviatoric_strain_rate

    @property
    def enr_shear_modulus(self):
        """
        Accessor on the right part of cracked cell shear modulus field
        """
        return self._enr_shear_modulus

    @property
    def enr_yield_stress(self):
        """
        Accessor on the right part of cracked cell yield stress field
        """
        return self._enr_yield_stress

    @property
    def enr_equivalent_plastic_strain_rate(self):
        """
        Accessor on the right part of cracked cell equivalent plastic strain rate at time t
        """
        return self._enr_equivalent_plastic_strain_rate

    @property
    def enr_plastic_strain_rate(self):
        """
        Accessor on the right part of cracked cell plastic strain rate tensor at time t
        """
        return self._enr_plastic_strain_rate

    @property
    def pressure_field(self):
        """
        :return: pressure field
        :rtype: np.array
        """
        return self.reconstruct_enriched_hydro_field(self.pressure, "enr_pressure")

    @property
    def density_field(self):
        """
        :return: density field
        :rtype: np.array
        """
        return self.reconstruct_enriched_hydro_field(self.density, "enr_density")

    @property
    def energy_field(self):
        """
        :return: energy field
        :rtype: np.array
        """
        return self.reconstruct_enriched_hydro_field(self.energy, "enr_energy")

    @property
    def artificial_viscosity_field(self):
        """
        :return: artificial viscosity field
        :rtype: np.array
        """
        return self.reconstruct_enriched_hydro_field(self.pseudo,
                                                     "enr_artificial_viscosity")

    @property
    def coordinates_field(self):
        """
        Returns the coordinates of cells, including left and right parts
        :return:
        """
        insertion_field = self._enr_coordinates_x[self.enriched]  # filtre(taille nb_disc)
        res = np.copy(self._coordinates_x[:, 0])
        offset = 0
        for enr_index in np.where(self.enriched)[0]:
            res = np.insert(res, enr_index + offset + 1, insertion_field[offset, 1])
            offset += 1
        return res

    @property
    def enr_coordinates_x(self):
        """
        Accesseur sur le tableau des cooordonnées enrichies
        """
        return self._enr_coordinates_x

    @property
    def stress_xx_field(self):
        """
        :return: sigma_xx field
        :rtype: np.array
        """
        return self.reconstruct_enriched_elasto_field(self.stress, "enr_stress")

    @property
    def deviatoric_stress_field(self):
        """
        :return: (deviateur de sigma)_xx field
        :rtype: np.array
        """
        return self.reconstruct_enriched_elasto_field(self.deviatoric_stress_current,
                                                      "enr_deviatoric_stress_current")

    @property
    def classical(self):
        """
        :return: a mask where True indicate a classical cell
        """
        return self._classical

    @property
    def enriched(self):
        """
        :return: a mask where True indicates an enrich cell
        """
        return ~self.classical

    def __str__(self):
        message = "<--ENRICHED CELLS COLLECTION-->" + os.linesep
        message += "Classical elements are:"
        message += str(self.classical) + os.linesep
        message += "Enriched elements are:"
        message += str(self.enriched)
        return message

    def print_infos(self):
        """
        Printing informations about Elements
            A REECRIRE AU PROPRE; NOTATIONS ONT CHANGE
        """
        message = "{}\n".format(self.__class__)
        for disc in Discontinuity.discontinuity_list():
            cell_i = disc.get_ruptured_cell_id
            message += "---- Discontinuity {:} ----".format(disc.label)
            # Density
            message += "==> masse volumique classique à t = {}\n". \
                        format(self.density.current_value[cell_i])
            message += "==> masse volumique enrichie à t = {}\n". \
                        format(self.enr_density.current_value[cell_i])
            message += "==> masse volumique classique à t+dt = {}\n". \
                        format(self.density.new_left_value[cell_i])
            message += "==> masse volumique enrichie à t+dt = {}\n". \
                        format(self.enr_density.new_value[cell_i])
            # Size of each part of the cracked cell
            message += "==> taille à gauche à t = {}\n". \
                format(disc.left_size.current_value)
            message += "==> taille à droite à t = {}\n". \
                format(disc.right_size.current_value)
            message += "==> taille à gauche à t+dt = {}\n". \
                format(disc.left_size.new_value)
            message += "==> taille à droite à t+dt = {}\n". \
                format(disc.right_size.new_value)
            # Pressure
            message += "==> pression à gauche à t = {}\n". \
                format(self.pressure.current_value[cell_i])
            message += "==> pression à droite à t = {}\n". \
                format(self.enr_pressure.current_value[cell_i])
            message += "==> pression à gauche à t+dt = {}\n". \
                format(self.pressure.new_value[cell_i])
            message += "==> pression à droite à t+dt = {}\n". \
                format(self.enr_pressure.new_value[cell_i])
            # Sound speed
            message += "==> vitesse du son à gauche à t = {}\n". \
                format(self.sound_velocity.current_value[cell_i])
            message += "==> vitesse du son à droite à t = {}\n". \
                format(self.enr_sound_velocity.current_value[cell_i])
            message += "==> vitesse du son à gauche à t+dt = {}\n". \
                format(self.sound_velocity.new_value[cell_i])
            message += "==> vitesse du son à droite à t+dt = {}\n". \
                format(self.enr_sound_velocity.new_value[cell_i])
            # Energy
            message += "==> énergie à gauche à t = {}\n". \
                format(self.energy.current_value[cell_i])
            message += "==> énergie à droite à t = {}\n". \
                format(self.enr_energy.current_value[cell_i])
            message += "==> énergie à gauche à t+dt = {}\n". \
                format(self.energy.new_value[cell_i])
            message += "==> énergie à droite à t+dt = {}\n". \
                format(self.enr_energy.new_value[cell_i])
            # Pseudo viscosity
            message += "==> pseudo à gauche = {}\n". \
                format(self.pseudo.current_value[cell_i])
            message += "==> pseudo à droite = {}\n". \
                format(self.enr_artificial_viscosity.current_value[cell_i])
        print(message)

    def compute_enriched_elements_new_pressure(self, delta_t):
        """
        Compute pressure, internal energy and sound velocity in left and right parts of
        the enriched elements

        :param delta_t: time step
        """
        target_model = self.data.material_target.constitutive_model
        # Fracture cannot occur on the projectile => check only the  target model to know if
        # elasticity or plasticity is activated
        elasticity_activated = (target_model.elasticity_model is not None)
        plasticity_activated = (target_model.plasticity_model is not None)

        mask = self.enriched
        if elasticity_activated or plasticity_activated:
            self.enr_energy.current_value[mask] += \
                OneDimensionCell.add_elastic_energy_method(
                    delta_t, self.enr_density.current_value[mask],
                    self.enr_density.new_value[mask],
                    self.enr_deviatoric_stress_current[mask],
                    self.enr_deviatoric_stress_new[mask],
                    self.enr_deviatoric_strain_rate[mask])

        # Initialize local parameters :
        density_right = self.enr_density.current_value[mask]
        density_right_new = self.enr_density.new_value[mask]
        pressure_right = self.enr_pressure.current_value[mask]
        pressure_right_new = self.enr_pressure.new_value[mask]
        energy_right = self.enr_energy.current_value[mask]
        energy_right_new = self.enr_energy.new_value[mask]
        pseudo_right = self.enr_artificial_viscosity.current_value[mask]
        cson_right_new = self.enr_sound_velocity.new_value[mask]
        # Call EOS :
        energy_new_right_value, pressure_new_right_value, sound_velocity_new_right_value = \
            OneDimensionCell.apply_equation_of_state(
                self, self._target_eos,
                density_right, density_right_new, pressure_right,
                pressure_right_new, energy_right, energy_right_new,
                pseudo_right, cson_right_new)

        # Save results :
        self.enr_pressure.new_value[mask] = pressure_new_right_value
        self.enr_energy.new_value[mask] = energy_new_right_value
        self.enr_sound_velocity.new_value[mask] = sound_velocity_new_right_value

    def compute_enriched_elements_new_part_size(self, time_step, node_velocity):
        """
        Compute the new size of each ruptured element part (left size and right size)

        :param time_step: time step
        :param node_velocity: array, node velocities
        """
        for disc in Discontinuity.discontinuity_list():
            u_left, u_right = (
                OneDimensionHansboEnrichedCell.compute_discontinuity_borders_velocity(
                    disc, node_velocity))
            u_node_left = node_velocity[disc.in_nodes]
            u_node_right = node_velocity[disc.out_nodes]
            self.left_part_size.new_value[disc.get_ruptured_cell_id] = (
                self.left_part_size.current_value[disc.get_ruptured_cell_id]
                + (u_left - u_node_left) * time_step)
            self.right_part_size.new_value[disc.get_ruptured_cell_id] = (
                self.right_part_size.current_value[disc.get_ruptured_cell_id]
                + (u_node_right - u_right) * time_step)

    def compute_enriched_elements_new_density(self):
        """
        Compute the new densities for left and right parts of the ruptured element
        (from mass conservation equation)
        """
        mask = self.enriched
        density_left = self.density.current_value[mask]
        density_right = self.enr_density.current_value[mask]
        size_left_current = self.left_part_size.current_value[mask]
        size_left_new = self.left_part_size.new_value[mask]
        size_right_current = self.right_part_size.current_value[mask]
        size_right_new = self.right_part_size.new_value[mask]

        self.density.new_value[mask] = density_left * size_left_current / size_left_new
        self.enr_density.new_value[mask] = (density_right *
                                            size_right_current / size_right_new)

    def compute_enriched_elements_new_pseudo(self, delta_t):
        """
        Compute the new artificial viscosity of the enriched_cells

        :param delta_t: time_step
        """
        mask = self.enriched
        if not mask.any():
            return
        # Left part :
        density_left = self.density.current_value[mask]
        density_left_new = self.density.new_value[mask]
        sound_velocity_left = self.sound_velocity.current_value[mask]
        size_left = self.left_part_size.new_value[mask]
        pseudo_left = OneDimensionCell.compute_pseudo(delta_t, density_left, density_left_new,
                                                      size_left, sound_velocity_left,
                                                      self.data.numeric.a_pseudo,
                                                      self.data.numeric.b_pseudo)
        # Right part :
        density_right = self.enr_density.current_value[mask]
        density_right_new = self.enr_density.new_value[mask]
        sound_velocity_right = self.enr_sound_velocity.current_value[mask]
        size_right = self.right_part_size.new_value[mask]
        pseudo_right = OneDimensionCell.compute_pseudo(delta_t, density_right, density_right_new,
                                                       size_right, sound_velocity_right,
                                                       self.data.numeric.a_pseudo,
                                                       self.data.numeric.b_pseudo)
        self.pseudo.new_value[mask] = pseudo_left
        self.enr_artificial_viscosity.new_value[mask] = pseudo_right

    def compute_enriched_stress_tensor(self):
        """
        Compute the complete enriched stress tensor : sigma = -(p+q) I + S
        """
        mask = self.enriched
        # Right part
        self.enr_stress[mask, 0] = \
            self.enr_deviatoric_stress_new[mask, 0] - \
            (self.enr_pressure.new_value[mask] +
             self.enr_artificial_viscosity.new_value[mask])
        self.enr_stress[mask, 1] = \
            self.enr_deviatoric_stress_new[mask, 1] - \
            (self.enr_pressure.new_value[mask] +
             self.enr_artificial_viscosity.new_value[mask])
        self.enr_stress[mask, 2] = \
            self.enr_deviatoric_stress_new[mask, 2] - \
            (self.enr_pressure.new_value[mask] +
             self.enr_artificial_viscosity.new_value[mask])

    def compute_enriched_deviatoric_strain_rate(self, dt: float,  # pylint: disable=invalid-name
                                                node_coord_new: np.array,
                                                node_velocity_new: np.array) -> None:
        """
        Compute the deviatoric strain rate for enriched cells

        :param dt: time step
        :param node_coord_new: array, new nodes coordinates
        :param node_velocity_new: array, new nodes velocity
        """
        disc_list = Discontinuity.discontinuity_list()
        if not disc_list:
            return

        nodes_in = Discontinuity.in_nodes.flatten()
        nodes_out = Discontinuity.out_nodes.flatten()
        mask_cells_arr = Discontinuity.ruptured_cell_id.flatten()
        eps_arr = Discontinuity.discontinuity_position.flatten()
        u2g_arr = Discontinuity.enr_velocity_new[:, 0].flatten()
        u1d_arr = Discontinuity.enr_velocity_new[:, 1].flatten()
        u_noeuds_new_in_arr = node_velocity_new[nodes_in]
        u_noeuds_new_out_arr = node_velocity_new[nodes_out]
        x_noeuds_new_in_arr = node_coord_new[nodes_in]
        x_noeuds_new_out_arr = node_coord_new[nodes_out]

        xg_new_arr = np.concatenate((
            x_noeuds_new_in_arr,
            x_noeuds_new_in_arr + self.left_part_size.new_value[mask_cells_arr][np.newaxis].T),
            axis=1)
        xd_new_arr = np.concatenate((
            x_noeuds_new_out_arr - self.right_part_size.new_value[mask_cells_arr][np.newaxis].T,
            x_noeuds_new_out_arr), axis=1)

        u_discg_new_arr, u_discd_new_arr = self._compute_discontinuity_borders_velocity(
            eps_arr, u_noeuds_new_in_arr[:, 0], u1d_arr, u2g_arr, u_noeuds_new_out_arr[:, 0])

        ug_new_arr = np.concatenate((u_noeuds_new_in_arr, u_discg_new_arr[np.newaxis].T), axis=1)
        ud_new_arr = np.concatenate((u_discd_new_arr[np.newaxis].T, u_noeuds_new_out_arr), axis=1)

        deviator_left = OneDimensionCell.general_method_deviator_strain_rate(dt, xg_new_arr, ug_new_arr)
        deviator_right = OneDimensionCell.general_method_deviator_strain_rate(dt, xd_new_arr, ud_new_arr)

        self._deviatoric_strain_rate[mask_cells_arr] = deviator_left
        self._enr_deviatoric_strain_rate[mask_cells_arr] = deviator_right

    def compute_enriched_deviatoric_stress_tensor(self, node_coord_new, node_velocity_new,
                                                  delta_t):
        """
        Compute the deviatoric part of the stress tensor

        :param node_coord_new: array, new nodes coordinates
        :param node_velocity_new: array, new nodes velocity
        :param delta_t: float, time step
        """
        self.compute_enriched_deviatoric_strain_rate(delta_t, node_coord_new, node_velocity_new)
        # Compute rotation rate tensor : W = 0 en 1D

        mask = self.enriched
        enr_cells = np.where(mask)[0]
        if not mask.any():
            return
        
        # Left part
        G = self.shear_modulus.new_value[enr_cells]  # pylint: disable=invalid-name
        D = self._deviatoric_strain_rate[enr_cells, :]
        S = self._deviatoric_stress_current[enr_cells, :]
        self._deviatoric_stress_new[enr_cells, 0] = S[:, 0] + 2. * G * D[:, 0] * delta_t
        self._deviatoric_stress_new[enr_cells, 1] = S[:, 1] + 2. * G * D[:, 1] * delta_t
        self._deviatoric_stress_new[enr_cells, 2] = S[:, 2] + 2. * G * D[:, 2] * delta_t

        # Right part
        Gd = self.enr_shear_modulus.new_value[enr_cells]  # pylint: disable=invalid-name
        Dd = self.enr_deviatoric_strain_rate[enr_cells, :]
        Sd = self.enr_deviatoric_stress_current[enr_cells, :]
        self._enr_deviatoric_stress_new[enr_cells, 0] = Sd[:, 0] + 2. * Gd * Dd[:, 0] * delta_t
        self._enr_deviatoric_stress_new[enr_cells, 1] = Sd[:, 1] + 2. * Gd * Dd[:, 1] * delta_t
        self._enr_deviatoric_stress_new[enr_cells, 2] = Sd[:, 2] + 2. * Gd * Dd[:, 2] * delta_t


    def compute_enriched_shear_modulus(self, shear_modulus_model):
        """
        Compute the shear modulus for ruptured cell

        :param shear_modulus_model: model to compute the shear modulus
        """
        mask = self.enriched
        if not mask.any():
            return
        self.enr_shear_modulus.new_value[mask] = \
            shear_modulus_model.compute(self.enr_density.new_value[mask])

    def apply_plasticity_enr(self, mask_mesh, delta_t):
        """
        Apply plasticity treatment if criterion is activated :

        - compute yield stress
        - tests plasticity criterion
        - compute plastic strain rate for plastic cells
        """
        mask = np.logical_and(self.plastic_enr_cells, mask_mesh)
        if not mask.any():
            return
        # Right part : right part of enriched cells is plastic ? => self.plastic_enr_cells
        invariant_j2_el = np.sqrt(compute_second_invariant(self.enr_deviatoric_stress_new[mask]))
        yield_stress = self.enr_yield_stress.new_value[mask]
        shear_modulus = self.enr_shear_modulus.new_value[mask]
        radial_return = self._compute_radial_return(invariant_j2_el, yield_stress)
        dev_stress = self._enr_deviatoric_stress_new[mask]
        self._enr_plastic_strain_rate[mask] = \
            self._compute_plastic_strain_rate_tensor(radial_return, shear_modulus,
                                                     delta_t, dev_stress)
        self._enr_equivalent_plastic_strain_rate[mask] = \
            self._compute_equivalent_plastic_strain_rate(invariant_j2_el, shear_modulus,
                                                         yield_stress, delta_t)
        self._enr_deviatoric_stress_new[mask] *= radial_return[np.newaxis].T

    def compute_enriched_yield_stress(self, yield_stress_model):
        """
        Compute the yield stress for ruptured cells

        :param yield_stress_model: model to compute the yield stress
        """
        mask = self.enriched
        self.enr_yield_stress.new_value[mask] = \
            yield_stress_model.compute(self.enr_density.new_value[mask])

    def compute_enriched_elements_new_time_step(self):
        """
        Compute the new local time step.
        The calculation is equivalent to a remeshing time step and thus underestimates the
        time step for the enriched cells
        """
        cfl = self.data.numeric.cfl
        cfl_pseudo = self.data.numeric.cfl_pseudo
        mask = self.enriched

        # Left part
        density_left = self.density.current_value[mask]
        density_left_new = self.density.new_value[mask]
        sound_velocity_left_new = self.sound_velocity.new_value[mask]
        pseudo_left = self.pseudo.current_value[mask]
        pseudo_left_new = self.pseudo.new_value[mask]
        dt_g = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, density_left, density_left_new,
                                                  self.left_part_size.new_value[mask],
                                                  sound_velocity_left_new,
                                                  pseudo_left, pseudo_left_new)
        # Right part
        density_right = self.enr_density.current_value[mask]
        density_right_new = self.enr_density.new_value[mask]
        sound_velocity_right_new = self.enr_sound_velocity.new_value[mask]
        pseudo_right = self.enr_artificial_viscosity.current_value[mask]
        pseudo_right_new = self.enr_artificial_viscosity.new_value[mask]
        dt_d = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, density_right, density_right_new,
                                                  self.right_part_size.new_value[mask],
                                                  sound_velocity_right_new,
                                                  pseudo_right, pseudo_right_new)
        if mask.any():
            self._dt[mask] = np.min(np.array([dt_g, dt_d]), axis=0)

    def cell_enr_increment(self):
        """
        Increment the enriched cell variables
        """
        # Thermodynamics
        self._enr_density.increment_values()
        self._enr_pressure.increment_values()
        self._enr_energy.increment_values()
        self._enr_artificial_viscosity.increment_values()
        self._enr_sound_velocity.increment_values()
        # Kinematics
        self._left_part_size.increment_values()
        self._right_part_size.increment_values()
        # Elasticity
        self._enr_deviatoric_stress_current[:] = self._enr_deviatoric_stress_new[:]

    def compute_new_coordinates(self, topology, node_coord):
        """
        Compute the coordinates of the cell center
        :param topology : mesh connectivity
        :param node_coord: coordinates of the nodes
        :return:
        """
        mask = self.enriched
        self._coordinates_x = node_coord[:-1, 0] + self.size_t_plus_dt / 2.
        cell_coord_from_right = node_coord[1:, 0] - self.size_t_plus_dt / 2.
        self._coordinates_x[mask] += (- self.size_t_plus_dt[mask]
                                      + self._left_part_size.new_value[mask]) / 2.
        self._enr_coordinates_x[mask] = cell_coord_from_right[mask] + \
            (self.size_t_plus_dt[mask] - self._right_part_size.new_value[mask]) / 2.
