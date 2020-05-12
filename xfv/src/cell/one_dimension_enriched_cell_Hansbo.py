# -*- coding: utf-8 -*-
"""
Implementing the Element1dEnriched class for Hansbo&Hansbo enrichment
"""
import os
import numpy as np

from xfv.src.cell.one_dimension_cell import OneDimensionCell
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.utilities.stress_invariants_calculation import compute_second_invariant
from xfv.src.fields.field import Field


# noinspection PyArgumentList
class OneDimensionHansboEnrichedCell(OneDimensionCell):  # pylint: disable=too-many-public-methods
    """
    A collection of 1d enriched elements. Treatment for Hansbo enrichment
    """
    @classmethod
    def compute_discontinuity_borders_velocity(cls, disc, node_velocity):
        """
        Compute the velocities of points at the discontinuity border
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
        ug = u2g * epsilon + u1g * (1. - epsilon)  # pylint: disable=invalid-name
        ud = u2d * epsilon + u1d * (1. - epsilon)  # pylint: disable=invalid-name
        return ug, ud

    def __init__(self, number_of_elements: int):
        """
        Build the class OneDimensionHansboEnrichedCell
        :param number_of_elements: total number of cells
        """
        super(OneDimensionHansboEnrichedCell, self).__init__(number_of_elements)
        #
        print(self._fields_manager)
        self._classical = np.ones(self.number_of_cells, dtype=np.bool, order='C')

        # Cell parts geometry
        self._left_part_size = Field(self.number_of_cells,
                                     current_value=np.zeros([self.number_of_cells]),
                                     new_value=np.zeros([self.number_of_cells]))
        self._right_part_size = Field(self.number_of_cells,
                                      current_value=np.zeros([self.number_of_cells]),
                                      new_value=np.zeros([self.number_of_cells]))

        # Additional dof: thermodynamics
        self._additional_dof_pressure = Field(self.number_of_cells,
                                              current_value=self.number_of_cells,
                                              new_value=np.zeros([self.number_of_cells]))
        self._additional_dof_density = Field(self.number_of_cells,
                                             current_value=np.zeros([self.number_of_cells]),
                                             new_value=np.zeros([self.number_of_cells]))
        self._additional_dof_energy = Field(self.number_of_cells,
                                            current_value=np.zeros([self.number_of_cells]),
                                            new_value=np.zeros([self.number_of_cells]))
        self._additional_dof_artificial_viscosity = \
            Field(self.number_of_cells, current_value=np.zeros([self.number_of_cells]),
                  new_value=np.zeros([self.number_of_cells]))
        self._additional_dof_sound_velocity = Field(self.number_of_cells,
                                                    current_value=np.zeros([self.number_of_cells]),
                                                    new_value=np.zeros([self.number_of_cells]))

        # Additional dof : elasticity
        self._additional_dof_stress = np.zeros([self.number_of_cells, 3])
        self._additional_dof_deviatoric_stress_current = np.zeros([self.number_of_cells, 3])
        self._additional_dof_deviatoric_stress_new = np.zeros([self.number_of_cells, 3])
        self._additional_dof_deviatoric_strain_rate = np.zeros([self.number_of_cells, 3])
        self._additional_dof_shear_modulus = Field(self.number_of_cells,
                                                   current_value=np.zeros([self.number_of_cells]),
                                                   new_value=np.zeros([self.number_of_cells]))
        self._additional_dof_yield_stress = Field(self.number_of_cells,
                                                  current_value=np.zeros([self.number_of_cells]),
                                                  new_value=np.zeros([self.number_of_cells]))
        self._additional_dof_equivalent_plastic_strain_rate = np.zeros([self.number_of_cells])
        self._additional_dof_plastic_strain_rate = np.zeros([self.number_of_cells, 3])
        self.plastic_enr_cells = np.zeros([self.number_of_cells],
                                          dtype=bool)  # Indicator right part of enr cell is plastic

    def initialize_additional_cell_dof(self, disc: Discontinuity):
        """
        Values to initialize the right part fields when discontinuity disc is created
        :param disc : the current discontinuity
        """
        mask = disc.mask_ruptured_cell
        # Initialization of the current field value
        self.additional_dof_density.current_value[mask] = \
            np.copy(self.density.current_value[mask])
        self.additional_dof_pressure.current_value[mask] = \
            np.copy(self.pressure.current_value[mask])
        self.additional_dof_sound_velocity.current_value[mask] = \
            np.copy(self.sound_velocity.current_value[mask])
        self.additional_dof_energy.current_value[mask] = np.copy(self.energy.current_value[mask])
        self.additional_dof_artificial_viscosity.current_value[mask] = \
            np.copy(self.pseudo.current_value[mask])
        self._additional_dof_deviatoric_stress_current[mask] = \
            np.copy(self._deviatoric_stress_current[mask])
        self.additional_dof_shear_modulus.current_value[mask] = \
            np.copy(self.shear_modulus.current_value[mask])
        self.additional_dof_yield_stress.current_value[mask] = \
            np.copy(self.yield_stress.current_value[mask])

        # Initialization of new value field
        # (so that the current value is not erased if the field is not updated in current step)
        self.additional_dof_density.new_value[mask] = np.copy(self.density.new_value[mask])
        self.additional_dof_pressure.new_value[mask] = np.copy(self.pressure.new_value[mask])
        self.additional_dof_sound_velocity.new_value[mask] = \
            np.copy(self.sound_velocity.new_value[mask])
        self.additional_dof_energy.new_value[mask] = np.copy(self.energy.new_value[mask])
        self.additional_dof_artificial_viscosity.new_value[mask] = \
            np.copy(self.pseudo.new_value[mask])
        self.additional_dof_shear_modulus.new_value[mask] = \
            np.copy(self.shear_modulus.new_value[mask])
        self.additional_dof_yield_stress.new_value[mask] = \
            np.copy(self.yield_stress.new_value[mask])
        self._additional_dof_deviatoric_stress_new[mask] = \
            np.copy(self._deviatoric_stress_new[mask])
        # Other quantities initialization
        self._additional_dof_deviatoric_strain_rate[mask] = \
            np.copy(self._deviatoric_strain_rate[mask])
        self._additional_dof_stress[mask] = np.copy(self._stress[mask])
        self._additional_dof_equivalent_plastic_strain_rate[mask] = \
            np.copy(self._equivalent_plastic_strain_rate[mask])

    def reconstruct_enriched_hydro_field(self, classical_field: Field, enriched_field_name: str):
        """
        True field reconstruction from the classical and enriched fields
        :param classical_field: classical field
        :param enriched_field_name: name of the enriched field
        :return: complete field
        :rtype np.array
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
        :rtype np.array
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
    def additional_dof_pressure(self):
        """
        Accessor on the additional cell pressure field
        """
        return self._additional_dof_pressure

    @property
    def additional_dof_density(self):
        """
        Accessor on the additional cell density field
        """
        return self._additional_dof_density

    @property
    def additional_dof_sound_velocity(self):
        """
        Accessor on the additional cell sound speed field
        """
        return self._additional_dof_sound_velocity

    @property
    def additional_dof_energy(self):
        """
        Accessor on the additional cell internal energy field
        """
        return self._additional_dof_energy

    @property
    def additional_dof_artificial_viscosity(self):
        """
        Accessor on the additional cell artificial viscosity field
        """
        return self._additional_dof_artificial_viscosity

    @property
    def additional_dof_stress(self):
        """
        Accessor on the additional cell stress at time t
        """
        return self._additional_dof_stress

    @property
    def additional_dof_stress_xx(self):
        """
        Accessor on the additional cell stress at time t
        """
        return self._additional_dof_stress[:, 0]

    @property
    def additional_dof_deviatoric_stress_current(self):
        """
        Accessor on the additional cell deviatoric stress at time t
        """
        return self._additional_dof_deviatoric_stress_current

    @property
    def additional_dof_deviatoric_stress_new(self):
        """
        Accessor on the additional cell deviatoric stress at time t+dt
        """
        return self._additional_dof_deviatoric_stress_new

    @property
    def additional_dof_deviatoric_strain_rate(self):
        """
        Accessor on the additional cell deviatoric strain rate at time t
        """
        return self._additional_dof_deviatoric_strain_rate

    @property
    def additional_dof_shear_modulus(self):
        """
        Accessor on the additional cell shear modulus field
        """
        return self._additional_dof_shear_modulus

    @property
    def additional_dof_yield_stress(self):
        """
        Accessor on the additional cell yield stress field
        """
        return self._additional_dof_yield_stress

    @property
    def additional_dof_equivalent_plastic_strain_rate(self):
        """
        Accessor on the additional cell equivalent plastic strain rate at time t
        """
        return self._additional_dof_equivalent_plastic_strain_rate

    @property
    def additional_dof_plastic_strain_rate(self):
        """
        Accessor on the additional cell plastic strain rate tensor at time t
        """
        return self._additional_dof_plastic_strain_rate

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
            cell_i = disc.ruptured_cell_id
            message += "---- Discontinuity {:} ----".format(disc.label)
            # Density
            message += "==> masse volumique classique à t = {}\n". \
                        format(self.density.current_value[cell_i])
            message += "==> masse volumique enrichie à t = {}\n". \
                        format(self.additional_dof_density.current_value[cell_i])
            message += "==> masse volumique classique à t+dt = {}\n". \
                        format(self.density.new_left_value[cell_i])
            message += "==> masse volumique enrichie à t+dt = {}\n". \
                        format(self.additional_dof_density.new_value[cell_i])
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
                format(self.additional_dof_pressure.current_value[cell_i])
            message += "==> pression à gauche à t+dt = {}\n". \
                format(self.pressure.new_value[cell_i])
            message += "==> pression à droite à t+dt = {}\n". \
                format(self.additional_dof_pressure.new_value[cell_i])
            # Sound speed
            message += "==> vitesse du son à gauche à t = {}\n". \
                format(self.sound_velocity.current_value[cell_i])
            message += "==> vitesse du son à droite à t = {}\n". \
                format(self.additional_dof_sound_velocity.current_value[cell_i])
            message += "==> vitesse du son à gauche à t+dt = {}\n". \
                format(self.sound_velocity.new_value[cell_i])
            message += "==> vitesse du son à droite à t+dt = {}\n". \
                format(self.additional_dof_sound_velocity.new_value[cell_i])
            # Energy
            message += "==> énergie à gauche à t = {}\n". \
                format(self.energy.current_value[cell_i])
            message += "==> énergie à droite à t = {}\n". \
                format(self.additional_dof_energy.current_value[cell_i])
            message += "==> énergie à gauche à t+dt = {}\n". \
                format(self.energy.new_value[cell_i])
            message += "==> énergie à droite à t+dt = {}\n". \
                format(self.additional_dof_energy.new_value[cell_i])
            # Pseudo viscosity
            message += "==> pseudo à gauche = {}\n". \
                format(self.pseudo.current_value[cell_i])
            message += "==> pseudo à droite = {}\n". \
                format(self.additional_dof_artificial_viscosity.current_value[cell_i])
        print(message)

    def compute_enriched_elements_new_pressure(self, delta_t):
        """
        Compute pressure, internal energy and sound velocity in left and right parts of
        the enriched elements
        :param delta_t : time step
        """
        target_model = self.data.material_target.constitutive_model
        # Fracture cannot occur on the projectile => check only the  target model to know if
        # elasticity or plasticity is activated
        elasticity_activated = (target_model.elasticity_model is not None)
        plasticity_activated = (target_model.plasticity_model is not None)

        mask = self.enriched
        if elasticity_activated or plasticity_activated:
            self.additional_dof_energy.current_value[mask] += \
                OneDimensionCell.add_elastic_energy_method(
                    delta_t, self.additional_dof_density.current_value[mask],
                    self.additional_dof_density.new_value[mask],
                    self.additional_dof_deviatoric_stress_current[mask],
                    self.additional_dof_deviatoric_stress_new[mask],
                    self.additional_dof_deviatoric_strain_rate[mask])

        # Initialize local parameters :
        density_right = self.additional_dof_density.current_value[mask]
        density_right_new = self.additional_dof_density.new_value[mask]
        pressure_right = self.additional_dof_pressure.current_value[mask]
        pressure_right_new = self.additional_dof_pressure.new_value[mask]
        energy_right = self.additional_dof_energy.current_value[mask]
        energy_right_new = self.additional_dof_energy.new_value[mask]
        pseudo_right = self.additional_dof_artificial_viscosity.current_value[mask]
        cson_right_new = self.additional_dof_sound_velocity.new_value[mask]
        # Call EOS :
        energy_new_right_value, pressure_new_right_value, sound_velocity_new_right_value = \
            OneDimensionCell.apply_equation_of_state(
                self, self._target_eos,
                density_right, density_right_new, pressure_right,
                pressure_right_new, energy_right, energy_right_new,
                pseudo_right, cson_right_new)

        # Save results :
        self.additional_dof_pressure.new_value[mask] = pressure_new_right_value
        self.additional_dof_energy.new_value[mask] = energy_new_right_value
        self.additional_dof_sound_velocity.new_value[mask] = sound_velocity_new_right_value

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
            self.left_part_size.new_value[disc.ruptured_cell_id] = \
                self.left_part_size.current_value[disc.ruptured_cell_id] + (ug - u1g) * time_step
            self.right_part_size.new_value[disc.ruptured_cell_id] = \
                self.right_part_size.current_value[disc.ruptured_cell_id] + (u2d - ud) * time_step

    def compute_enriched_elements_new_density(self):
        """
        Compute the new densities for left and right parts of the ruptured element
        (from mass conservation equation)
        """
        mask = self.enriched
        density_left = self.density.current_value[mask]
        density_right = self.additional_dof_density.current_value[mask]
        size_left_current = self.left_part_size.current_value[mask]
        size_left_new = self.left_part_size.new_value[mask]
        size_right_current = self.right_part_size.current_value[mask]
        size_right_new = self.right_part_size.new_value[mask]

        self.density.new_value[mask] = density_left * size_left_current / size_left_new
        self.additional_dof_density.new_value[mask] = (density_right *
                                                       size_right_current / size_right_new)

    def compute_enriched_elements_new_pseudo(self, delta_t):
        """
        Compute the new artificial viscosity of the enriched_cells
        :param delta_t: time_step
        """
        mask = self.enriched
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
        density_right = self.additional_dof_density.current_value[mask]
        density_right_new = self.additional_dof_density.new_value[mask]
        sound_velocity_right = self.additional_dof_sound_velocity.current_value[mask]
        size_right = self.right_part_size.new_value[mask]
        pseudo_right = OneDimensionCell.compute_pseudo(delta_t, density_right, density_right_new,
                                                       size_right, sound_velocity_right,
                                                       self.data.numeric.a_pseudo,
                                                       self.data.numeric.b_pseudo)
        self.pseudo.new_value[mask] = pseudo_left
        self.additional_dof_artificial_viscosity.new_value[mask] = pseudo_right

    def compute_enriched_stress_tensor(self):
        """
        Compute the complete enriched stress tensor : sigma = -(p+q) I + S
        """
        mask = self.enriched
        # Right part
        self.additional_dof_stress[mask, 0] = \
            self.additional_dof_deviatoric_stress_new[mask, 0] - \
            (self.additional_dof_pressure.new_value[mask] +
             self.additional_dof_artificial_viscosity.new_value[mask])
        self.additional_dof_stress[mask, 1] = \
            self.additional_dof_deviatoric_stress_new[mask, 1] - \
            (self.additional_dof_pressure.new_value[mask] +
             self.additional_dof_artificial_viscosity.new_value[mask])
        self.additional_dof_stress[mask, 2] = \
            self.additional_dof_deviatoric_stress_new[mask, 2] - \
            (self.additional_dof_pressure.new_value[mask] +
             self.additional_dof_artificial_viscosity.new_value[mask])

    def compute_enriched_deviatoric_strain_rate(self, dt,  # pylint: disable=invalid-name
                                                node_coord_new,
                                                node_velocity_new):
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
                               x_noeuds_new[0] + self.left_part_size.new_value[mask_cells]])
            xg_new = xg_new.reshape(1, 2)
            ug_new = np.array([u_noeuds_new[0], u_discg_new]).reshape(1, 2)
            # Compute the deviatoric strain rate tensor for left part
            D = OneDimensionCell.general_method_deviator_strain_rate(dt, xg_new, ug_new)  # np.array(True) to be consistent
            self._deviatoric_strain_rate[mask_cells] = D 

            # Creation of structure left - right data to call general_method_deviator_strain_rate
            # Left part cell : right boundary of discontinuity - node_right
            xd_new = np.array([x_noeuds_new[1] - self.right_part_size.new_value[mask_cells],
                               x_noeuds_new[1]])
            xd_new = xd_new.reshape(1, 2)
            ud_new = np.array([u_discd_new, u_noeuds_new[1]]).reshape(1, 2)
            D = OneDimensionCell.general_method_deviator_strain_rate(dt, xd_new, ud_new)
            # Compute the deviatoric strain rate tensor for right part
            self._additional_dof_deviatoric_strain_rate[mask_cells] = D

    def compute_enriched_deviatoric_stress_tensor(self, node_coord_new, node_velocity_new,
                                                  dt):  # pylint: disable=invalid-name
        """
        Compute the deviatoric part of the stress tensor
        :param node_coord_new : array, new nodes coordinates
        :param node_velocity_new : array, new nodes velocity
        :param dt : float, time step
        """
        self.compute_enriched_deviatoric_strain_rate(dt, node_coord_new, node_velocity_new)
        # Compute rotation rate tensor : W = 0 en 1D

        # Left part
        mask = self.enriched
        if not mask.any():
            return
        G = self.shear_modulus.new_value[mask]  # pylint: disable=invalid-name
        self._deviatoric_stress_new[mask, 0] = self._deviatoric_stress_current[mask, 0] + \
            2. * G * self._deviatoric_strain_rate[mask, 0] * dt
        self._deviatoric_stress_new[mask, 1] = self._deviatoric_stress_current[mask, 1] + \
            2. * G * self._deviatoric_strain_rate[mask, 1] * dt
        self._deviatoric_stress_new[mask, 2] = self._deviatoric_stress_current[mask, 2] + \
            2. * G * self._deviatoric_strain_rate[mask, 2] * dt
        # Right part
        Gd = self.additional_dof_shear_modulus.new_value[mask]  # pylint: disable=invalid-name
        self._additional_dof_deviatoric_stress_new[mask, 0] = \
            self.additional_dof_deviatoric_stress_current[mask, 0] + \
            2. * Gd * self.additional_dof_deviatoric_strain_rate[mask, 0] * dt
        self._additional_dof_deviatoric_stress_new[mask, 1] = \
            self.additional_dof_deviatoric_stress_current[mask, 1] + \
            2. * Gd * self.additional_dof_deviatoric_strain_rate[mask, 1] * dt
        self._additional_dof_deviatoric_stress_new[mask, 2] = \
            self.additional_dof_deviatoric_stress_current[mask, 2] +\
            2. * Gd * self.additional_dof_deviatoric_strain_rate[mask, 2] * dt

    def compute_enriched_shear_modulus(self, shear_modulus_model):
        """
        Compute the shear modulus for ruptured cell
        :param shear_modulus_model : model to compute the shear modulus
        """
        mask = self.enriched
        if not mask.any():
            return
        self.additional_dof_shear_modulus.new_value[mask] = \
            shear_modulus_model.compute(self.additional_dof_density.new_value[mask])

    def apply_plastic_correction_on_enriched_deviatoric_stress_tensor(self, mask_mesh):
        """
        Correct the elastic trial of deviatoric stress tensor when plasticity criterion is activated
        :param mask_mesh:  mask to identify the part of the mesh (projectile or target)
        """
        mask = np.logical_and(self.plastic_enr_cells, mask_mesh)
        # Right part of the cracked cell :
        invariant_j2_el_right = np.sqrt(compute_second_invariant(
            self.additional_dof_deviatoric_stress_new[mask]))
        yield_stress = self.additional_dof_yield_stress.new_value[mask]
        radial_return_right = yield_stress / invariant_j2_el_right
        self._additional_dof_deviatoric_stress_new[mask, 0] *= radial_return_right
        self._additional_dof_deviatoric_stress_new[mask, 1] *= radial_return_right
        self._additional_dof_deviatoric_stress_new[mask, 2] *= radial_return_right

    def compute_enriched_equivalent_plastic_strain_rate(self, mask_mesh,
                                                        dt):  # pylint: disable=invalid-name
        """
        Compute the plastic strain rate
        :param dt : time step
        :param mask_mesh:  mask to identify the part of the mesh (projectile or target)
        """
        mask = np.logical_and(self.plastic_enr_cells, mask_mesh)
        # Right part :
        invariant_j2_el_right = np.sqrt(compute_second_invariant(
            self.additional_dof_deviatoric_stress_new[mask]))  # elastic predictor
        shear_mod_right = self.additional_dof_shear_modulus.new_value[mask]
        yield_stress_right = self.additional_dof_yield_stress.new_value[mask]
        self._additional_dof_equivalent_plastic_strain_rate[mask] = \
            (invariant_j2_el_right - yield_stress_right) / (3. * shear_mod_right * dt)

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
        invariant_j2_el = np.sqrt(compute_second_invariant(self.additional_dof_deviatoric_stress_new[mask]))
        yield_stress = self.additional_dof_yield_stress.new_value[mask]
        shear_modulus = self.additional_dof_shear_modulus.new_value[mask]
        radial_return = self._compute_radial_return(invariant_j2_el, yield_stress)
        dev_stress = self.additional_dof_deviatoric_stress_new[mask]
        self._plastic_strain_rate[mask] = self._compute_plastic_strain_rate_tensor(radial_return, shear_modulus, delta_t, dev_stress)
        self._equivalent_plastic_strain_rate[mask] = self._compute_equivalent_plastic_strain_rate(invariant_j2_el, shear_modulus, yield_stress, delta_t)
        self._deviatoric_stress_new[mask] *= radial_return[np.newaxis].T


    def compute_enriched_plastic_strain_rate(self, mask_mesh, dt):  # pylint: disable=invalid-name
        """
        Compute the plastic strain rate tensor from elastic prediction and radial return
        (normal law for Von Mises plasticity) in cracked cells
        :param mask_mesh : mask to identify the part of the mesh (projectile or target)
        :param dt: time step
        """
        mask = np.logical_and(self.plastic_enr_cells, mask_mesh)
        # Right part : right part of enriched cells is plastic ? => self.plastic_enr_cells
        invariant_j2_el_right = np.sqrt(
            compute_second_invariant(self.additional_dof_deviatoric_stress_new[mask]))
        shear_mod_right = self.additional_dof_shear_modulus.new_value[mask]
        yield_stress_right = self.additional_dof_yield_stress.new_value[mask]
        radial_return_right = yield_stress_right / invariant_j2_el_right
        for i in range(0, 3):
            self._additional_dof_plastic_strain_rate[mask, i] = \
                (1 - radial_return_right) / (radial_return_right * 3 * shear_mod_right * dt) * \
                self._additional_dof_deviatoric_stress_new[mask, i]

    def compute_enriched_yield_stress(self, yield_stress_model):
        """
        Compute the yield stress for ruptured cells
        :param yield_stress_model : model to compute the yield stress
        """
        mask = self.enriched
        self.additional_dof_yield_stress.new_value[mask] = \
            yield_stress_model.compute(self.additional_dof_density.new_value[mask])

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
        density_right = self.additional_dof_density.current_value[mask]
        density_right_new = self.additional_dof_density.new_value[mask]
        sound_velocity_right_new = self.additional_dof_sound_velocity.new_value[mask]
        pseudo_right = self.additional_dof_artificial_viscosity.current_value[mask]
        pseudo_right_new = self.additional_dof_artificial_viscosity.new_value[mask]
        dt_d = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, density_right, density_right_new,
                                                  self.right_part_size.new_value[mask],
                                                  sound_velocity_right_new,
                                                  pseudo_right, pseudo_right_new)
        if mask.any():
            self._dt[mask] = np.min(np.array([dt_g, dt_d]), axis=0)

    def cell_additional_dof_increment(self):
        """
        Increment the enriched cell variables
        """
        # Thermodynamics
        self._additional_dof_density.increment_values()
        self._additional_dof_pressure.increment_values()
        self._additional_dof_energy.increment_values()
        self._additional_dof_artificial_viscosity.increment_values()
        self._additional_dof_sound_velocity.increment_values()
        # Kinematics
        self._left_part_size.increment_values()
        self._right_part_size.increment_values()
        # Elasticity
        self._additional_dof_deviatoric_stress_current[:] = \
            self._additional_dof_deviatoric_stress_new[:]
