#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A module implementing the Discontinuity class
"""
import numpy as np
from xfv.src.data.data_container import DataContainer
from xfv.src.fields.field import Field
from xfv.src.mass_matrix.one_dimension_enriched_mass_matrix_Hansbo import OneDimensionHansboEnrichedMassMatrix


class Discontinuity(object):
    """
    A class describing a discontinuity 1D
    """

    # A list of discontinuities
    __discontinuity_list = []

    def __init__(self, mask_in_nodes, mask_out_nodes, discontiuity_position_in_ruptured_element):
        """
        Initializing a single discontinuity after enrichment.
        :param mask_in_nodes: enriched nodes on the left of the discontinuity
        :param mask_out_nodes: enriched nodes on the right of the discontinuity
        :param discontiuity_position_in_ruptured_element ; position of the discontiuity in the ruptured element
        """
        # V�rification d'usage :
        if discontiuity_position_in_ruptured_element > 1 or discontiuity_position_in_ruptured_element < 0:
            raise ValueError("""Discontinuity position in cracked cell must be between 0 and 1""")

        if np.array([a in np.where(mask_out_nodes)[0] for a in np.where(mask_in_nodes)[0]]).any():
            raise ValueError("""A node cannot be both inside and outside the discontinuity""")

        Discontinuity.__discontinuity_list.append(self)
        self.__label = len(Discontinuity.__discontinuity_list)
        self.__mask_in_nodes = mask_in_nodes
        self.__mask_out_nodes = mask_out_nodes
        self.__mask_ruptured_cell = np.zeros(len(self.mask_in_nodes)-1, dtype=bool)
        self._ruptured_cell_id = 0
        self._discontinuity_position = discontiuity_position_in_ruptured_element
        self.__mass_matrix_updated = False
        self.__initialisation = False
        print("Building discontinuity number {:d}".format(self.__label))
        # Additional dof representing either the enriched Heaviside value or
        # the field value in the right part of enriched element.
        self._additional_dof_velocity_current = np.zeros([2, 1])
        self._additional_dof_velocity_new = np.zeros([2, 1])
        self._additional_dof_force = np.zeros([2, 1])

        # Additional dof: thermo
        self._left_part_size = Field(1, current_value=0., new_value=0.)
        self._right_part_size = Field(1, current_value=0., new_value=0.)
        self._additional_dof_pressure = Field(1, current_value=np.array([0.]), new_value=np.array([0.]))
        self._additional_dof_density = Field(1, current_value=np.array([0.]), new_value=np.array([0.]))
        self._additional_dof_energy = Field(1, current_value=np.array([0.]), new_value=np.array([0.]))
        self._additional_dof_artificial_viscosity = Field(1, current_value=np.array([0.]), new_value=np.array([0.]))
        self._additional_dof_sound_velocity = Field(1, current_value=np.array([0.]), new_value=np.array([0.]))

        # Additional dof : elasticity
        self._additional_dof_stress = np.zeros([1, 3])
        self._additional_dof_deviatoric_stress_current = np.zeros([1, 3])
        self._additional_dof_deviatoric_stress_new = np.zeros([1, 3])
        self._additional_dof_deviatoric_strain_rate = np.zeros([1, 3])
        self._additional_dof_shear_modulus = Field(1, current_value=np.array([0.]), new_value=np.array([0.]))
        self._additional_dof_yield_stress = Field(1, current_value=np.array([0.]), new_value=np.array([0.]))
        self._additional_dof_equivalent_plastic_strain_rate = 0.
        self._additional_dof_plastic_strain_rate = np.zeros([1, 3])
        self.plastic_cells = False  # La partie droite de la cell rompue est plastique

        # endommagement avec mod�le coh�sif (cr�� tout le temps mais reste � 0 si pas d'endommagement dans XDATA.json...
        target_dmg_model = DataContainer().material_target.damage_model
        if target_dmg_model and target_dmg_model.cohesive_model is not None:
            cohesive_strength = DataContainer().material_target.damage_model.cohesive_model.cohesive_strength
        else:
            cohesive_strength = 0
        self.cohesive_force = Field(1, current_value=cohesive_strength, new_value=cohesive_strength)
        self.discontinuity_opening = Field(1, current_value=0., new_value=0.)
        self.damage_variable = Field(1, current_value=0., new_value=0.)
        self.history_max_opening = 0.

        if target_dmg_model and target_dmg_model.cohesive_model is not None:
            self.history_min_cohesive_force = cohesive_strength
        else:  # dans le cas o�  pas de czm activ� : czm_model= None
            self.history_min_cohesive_force = 0.

        # cr�ation de la matrice de masse associ�e � la discontinuit�
        if DataContainer().material_target.failure_model.type_of_enrichment == "Hansbo":
            self.mass_matrix_enriched = OneDimensionHansboEnrichedMassMatrix(
                lump=DataContainer().material_target.failure_model.lump_mass_matrix)

    @classmethod
    def discontinuity_number(cls):
        return len(Discontinuity.__discontinuity_list)

    @classmethod
    def discontinuity_list(cls):
        return Discontinuity.__discontinuity_list

    @classmethod
    def get_discontinuity_associated_with_cell(cls, cell_id):
        """
        Parcourt lacollection des discontinuit�s � la recherche de celle qui contient la cell_id en argument
        :param cell_id: cell id � retrouver
        :return: Discontinuity qui matche
        """
        if cell_id == -1: return None
        try_index = 0
        while try_index < Discontinuity.discontinuity_number():
            disc = Discontinuity.discontinuity_list()[try_index]
            if disc._ruptured_cell_id == cell_id:
                return disc
            try_index += 1
        return None

    @property
    def position_in_ruptured_element(self):
        """
        Accessor on the position of the discontinuity in ruptured element
        """
        return self._discontinuity_position

    @property
    def label(self):
        """
        Accessor on label variable
        :return: label
        """
        return self.__label

    @property
    def mask_in_nodes(self):
        """
        Accessor on the mask on the nodes "in" the discontinuity
        :return: the mask on the nodes "in" the discontinuity
        """
        return self.__mask_in_nodes

    @property
    def mask_out_nodes(self):
        """
        Accessor on the mask on the nodes "out" the discontinuity
        :return: the mask on the nodes "out" the discontinuity
        """
        return self.__mask_out_nodes

    @property
    def mask_disc_nodes(self):
        """
        Accessor on the mask on the nodes of the discontinuity
        :return: the mask on the nodes "concerned" by the discontinuity
        """
        return np.logical_or(self.__mask_in_nodes, self.__mask_out_nodes)

    @property
    def mass_matrix_updated(self):
        """
        Accessor on the boolean that indicates if the mass matrix has been computed for this discontinuity

        :return: the boolean that indicates if the mass matrix has been computed for this discontinuity
        """
        return self.__mass_matrix_updated

    @property
    def initialisation(self):
        """
        Boolean to follow if the right part values have been correctly initialized for Hansbo enrichment.
        Indeed, those variables have physical meaning and are not zero
        Initialization is done with left part field values
        :return:
        """
        return self.__initialisation

    @property
    def ruptured_cell_id(self):
        """
        Returns the id of ruptured cell for the discontinuity
        """
        return int(self._ruptured_cell_id)

    @property
    def mask_ruptured_cell(self):
        """
        Returns a mask to identify the ruptured cell associated with discontinuity
        """
        return self.__mask_ruptured_cell

    def have_dof_been_initialized(self):
        """
        Set the __initialisation boolean to True
        """
        self.__initialisation = True

    def has_mass_matrix_been_computed(self):
        """
        Set the __mass_matrix_updated boolean to True
        """
        self.__mass_matrix_updated = True

    def find_ruptured_cell_id(self, topology):
        """
        Compute the cell id of ruptured element associated with this discontinuity
        self.ruptured_cell_id is an integer
        self.mask_ruptured_cell is an array of boolean = True for ruptured cell
        """
        self._ruptured_cell_id = topology.cells_in_contact_with_node[self.mask_in_nodes][0][1]
        self.__mask_ruptured_cell[self._ruptured_cell_id] = True


    def compute_discontinuity_new_opening(self, node_position):
        """

        :param node_position:
        :return:
        """
        xd_new = node_position[self.mask_out_nodes] - self.right_part_size.new_value
        xg_new = node_position[self.mask_in_nodes] + self.left_part_size.new_value
        self.discontinuity_opening.new_value = (xd_new - xg_new)[0][0]

    @property
    def discontinuity_position(self):
        """

        :return:
        """
        return self._discontinuity_position

    @property
    def additional_dof_velocity_current(self):
        return self._additional_dof_velocity_current

    @property
    def additional_dof_velocity_new(self):
        return self._additional_dof_velocity_new

    @property
    def additional_dof_force(self):
        return self._additional_dof_force

    @property
    def additional_dof_pressure(self):
        return self._additional_dof_pressure

    @property
    def additional_dof_density(self):
        return self._additional_dof_density

    @property
    def additional_dof_sound_velocity(self):
        return self._additional_dof_sound_velocity

    @property
    def additional_dof_energy(self):
        return self._additional_dof_energy

    @property
    def additional_dof_artificial_viscosity(self):
        return self._additional_dof_artificial_viscosity

    @property
    def additional_dof_stress(self):
        return self._additional_dof_stress

    @property
    def additional_dof_deviatoric_stress_current(self):
        return self._additional_dof_deviatoric_stress_current

    @property
    def additional_dof_deviatoric_stress_new(self):
        return self._additional_dof_deviatoric_stress_new

    @property
    def additional_dof_deviatoric_strain_rate(self):
        return self._additional_dof_deviatoric_strain_rate

    @property
    def additional_dof_shear_modulus(self):
        return self._additional_dof_shear_modulus

    @property
    def additional_dof_yield_stress(self):
        return self._additional_dof_yield_stress

    @property
    def additional_dof_equivalent_plastic_strain_rate(self):
        return self._additional_dof_equivalent_plastic_strain_rate

    @property
    def additional_dof_plastic_strain_rate(self):
        return self._additional_dof_plastic_strain_rate

    @property
    def left_part_size(self):
        return self._left_part_size

    @property
    def right_part_size(self):
        return self._right_part_size

    def additional_dof_increment(self):
        self._additional_dof_density.increment_values()
        self._additional_dof_pressure.increment_values()
        self._additional_dof_energy.increment_values()
        self._additional_dof_artificial_viscosity.increment_values()
        self._additional_dof_sound_velocity.increment_values()
        self._additional_dof_velocity_current[:] = self.additional_dof_velocity_new[:]
        self._left_part_size.increment_values()
        self._right_part_size.increment_values()
        self._additional_dof_deviatoric_stress_current[:] = self._additional_dof_deviatoric_stress_new[:]
        self.cohesive_force.increment_values()
        self.damage_variable.increment_values()
        self.discontinuity_opening.increment_values()


