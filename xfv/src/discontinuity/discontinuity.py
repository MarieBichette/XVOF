# -*- coding: utf-8 -*-
"""
A module implementing the Discontinuity class
"""
import numpy as np
from xfv.src.fields.field import Field
from xfv.src.data.enriched_mass_matrix_props import EnrichedMassMatrixProps


class Discontinuity:
    """
    A class describing a discontinuity 1D
    """

    # A list of discontinuities
    __discontinuity_list = []

    # Enriched variables
    enr_velocity_current = np.zeros([], dtype=float)
    enr_velocity_new = np.zeros([], dtype=float)
    enr_coordinates_current = np.zeros([], dtype=float)
    enr_coordinates_new = np.zeros([], dtype=float)
    enr_force = np.zeros([], dtype=float)

    # Information about the discontinuities
    discontinuity_position = np.zeros([], dtype=float)
    ruptured_cell_id = np.zeros([], dtype=int)
    in_nodes = np.zeros([], dtype=int)
    out_nodes = np.zeros([], dtype=int)

    def __init__(self, cell_id: int, mask_in_nodes: np.array, mask_out_nodes: np.array,
                 discontinuity_position_in_ruptured_element: float,
                 enriched_mass_matrix_props: EnrichedMassMatrixProps):
        """
        Initializing a single discontinuity after enrichment.

        :param cell_id: id of the cracked cell
        :param mask_in_nodes: enriched nodes on the left of the discontinuity
        :param mask_out_nodes: enriched nodes on the right of the discontinuity
        :param discontinuity_position_in_ruptured_element: position of the discontinuity in
        the ruptured element
        """
        # Some verification:
        if discontinuity_position_in_ruptured_element > 1 or \
                discontinuity_position_in_ruptured_element < 0:
            raise ValueError("""Discontinuity position in cracked cell must be between 0 and 1""")

        if np.array([a in np.where(mask_out_nodes)[0] for a in np.where(mask_in_nodes)[0]]).any():
            raise ValueError("""A node cannot be both inside and outside the discontinuity""")

        # Discontinuity registration
        Discontinuity.__discontinuity_list.append(self)
        self.__label = len(Discontinuity.__discontinuity_list)
        print("Building discontinuity number {:d}".format(self.__label))
        init = np.zeros((1, 2, 1))
        if self.__label == 1:
            Discontinuity.enr_velocity_current = np.copy(init)
            Discontinuity.enr_velocity_new = np.copy(init)
            Discontinuity.enr_coordinates_current = np.copy(init)
            Discontinuity.enr_coordinates_new = np.copy(init)
            Discontinuity.enr_force = np.copy(init)
            Discontinuity.discontinuity_position = np.zeros((1, 1))
            Discontinuity.ruptured_cell_id = np.zeros((1, 1), dtype=int)
            Discontinuity.in_nodes = np.zeros((1, 1), dtype=int)
            Discontinuity.out_nodes = np.zeros((1, 1), dtype=int)
            Discontinuity.critical_strength = np.zeros((1,1), dtype = float)
            Discontinuity.critical_separation = np.zeros((1,1), dtype = float)
        else:
            Discontinuity.enr_velocity_current = np.append(
                Discontinuity.enr_velocity_current, init, axis=0)
            Discontinuity.enr_velocity_new = np.append(
                Discontinuity.enr_velocity_new, init, axis=0)
            Discontinuity.enr_coordinates_current = np.append(
                Discontinuity.enr_coordinates_current, init, axis=0)
            Discontinuity.enr_coordinates_new = np.append(
                Discontinuity.enr_coordinates_new, init, axis=0)
            Discontinuity.enr_force = np.append(Discontinuity.enr_force, init, axis=0)
            Discontinuity.discontinuity_position = np.append(
                Discontinuity.discontinuity_position, np.zeros((1, 1)), axis=0)
            Discontinuity.ruptured_cell_id = np.append(
                Discontinuity.ruptured_cell_id, np.zeros((1, 1), dtype=int), axis=0)
            Discontinuity.in_nodes = np.append(
                Discontinuity.in_nodes, np.zeros((1, 1), dtype=int), axis=0)
            Discontinuity.out_nodes = np.append(
                Discontinuity.out_nodes, np.zeros((1, 1), dtype=int), axis=0)
            Discontinuity.critical_strength = np.append(
                Discontinuity.critical_strength, np.zeros((1, 1), dtype=float), axis=0)
            Discontinuity.critical_separation = np.append(
                Discontinuity.critical_separation, np.zeros((1, 1), dtype=float), axis=0)

        for ind, disc in enumerate(Discontinuity.discontinuity_list()):
            disc.enr_velocity_current = Discontinuity.enr_velocity_current[ind]
            disc.enr_velocity_new = Discontinuity.enr_velocity_new[ind]
            disc.enr_coordinates_current = Discontinuity.enr_coordinates_current[ind]
            disc.enr_coordinates_new = Discontinuity.enr_coordinates_new[ind]
            disc.enr_force = Discontinuity.enr_force[ind]
            disc.discontinuity_position = Discontinuity.discontinuity_position[ind]
            disc.ruptured_cell_id = Discontinuity.ruptured_cell_id[ind]
            disc.in_nodes = Discontinuity.in_nodes[ind]
            disc.out_nodes = Discontinuity.out_nodes[ind]
            disc.critical_strength = Discontinuity.critical_strength[ind]
            disc.critical_separation = Discontinuity.critical_separation[ind]
            


        self.__mask_in_nodes = mask_in_nodes
        self.in_nodes[:] = np.where(self.__mask_in_nodes)[0]
        self.__mask_out_nodes = mask_out_nodes
        self.out_nodes[:] = np.where(self.__mask_out_nodes)[0]

        # Save the cracked cell id associated to the disc
        # TODO : pourquoi array ??????
        self.ruptured_cell_id[:] = cell_id

        # Discontinuity cell information
        self.discontinuity_position[:] = discontinuity_position_in_ruptured_element

        # Indicators of discontinuity state
        self.__mass_matrix_updated = False

        # Damage indicators with cohesive zone model
        # (Always created but null if no damage data in the XDATA.json file...)
        self.cohesive_force = Field(1, current_value=0., new_value=0.)
        self.discontinuity_opening = Field(1, current_value=0., new_value=0.)
        self.dissipated_energy = Field(1, current_value=0., new_value=0.)
        self.damage_variable = Field(1, current_value=0., new_value=0.)
        self.history_max_opening = 0.
        self.history_min_cohesive_force = 1.e+30

        # Creation of the enriched mass matrix
        self.mass_matrix_enriched = enriched_mass_matrix_props.build_enriched_mass_matrix_obj()

    @classmethod
    def discontinuity_number(cls):
        """
        Returns the number of existing discontinuities
        """
        return len(Discontinuity.__discontinuity_list)

    @classmethod
    def discontinuity_list(cls):
        """
        Returns the list of all existing discontinuities
        """
        return Discontinuity.__discontinuity_list

    @classmethod
    def get_discontinuity_associated_with_cell(cls, cell_id: int):
        """
        Loop on the discontinuities collection to find the one that contain the cell cell_id

        :param cell_id: cell id to find
        :return: Discontinuity or None if no disc is found
        """
        if cell_id == -1:
            return None

        try_index = 0
        while try_index < Discontinuity.discontinuity_number():
            disc = Discontinuity.discontinuity_list()[try_index]
            if disc.get_ruptured_cell_id() == cell_id:
                return disc
            try_index += 1
        return None

    @property
    def position_in_ruptured_element(self):
        """
        Accessor on the relative position of the discontinuity in ruptured element
        """
        return self.discontinuity_position[0]

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
        Accessor on the boolean that indicates if the mass matrix has been computed

        :return: the boolean that indicates if the mass matrix has been computed
        """
        return self.__mass_matrix_updated

    @property
    def get_ruptured_cell_id(self):
        """
        Returns the id of ruptured cell for the discontinuity
        """
        return int(self.ruptured_cell_id[0])

    def has_mass_matrix_been_computed(self):
        """
        Set the __mass_matrix_updated boolean to True
        """
        self.__mass_matrix_updated = True

    def compute_discontinuity_new_opening(self, node_position: np.array):
        """
        Compute the discontinuity opening

        :param node_position: coordinates of the nodes
        """
        epsilon = self.discontinuity_position
        coord_g = node_position[self.mask_in_nodes]
        coord_d = node_position[self.mask_out_nodes]
        enr_coord_g = self.enr_coordinates_new[0]  # x2-
        enr_coord_d = self.enr_coordinates_new[1]  # x1+
        xg_new = (1 - epsilon) * coord_g + epsilon * enr_coord_g
        xd_new = (1 - epsilon) * enr_coord_d + epsilon * coord_d
        self.discontinuity_opening.new_value = (xd_new - xg_new)[0][0]

    def compute_critical_value(stress, energy):
        for ind in range(len(Discontinuity.critical_strength)):
            if abs(Discontinuity.critical_strength[ind]) < 1.e-16:
                Discontinuity.critical_strength[ind] = abs(stress[Discontinuity.ruptured_cell_id[ind]])
            if abs(Discontinuity.critical_separation[ind]) < 1.e-16:
                Discontinuity.critical_separation[ind] = 2*energy[Discontinuity.ruptured_cell_id[ind]]/Discontinuity.critical_strength[ind]
        
    def reinitialize_kinematics_after_contact(self):
        """
        Set the new velocity to the old one to cancel the increment that has lead to contact
        """
        self.enr_velocity_new[:] = self.enr_velocity_current[:]
        self.enr_coordinates_new[:] = self.enr_coordinates_current[:]

    def enr_increment(self):
        """
        Increment the variables of discontinuity
        """
        # Kinematics
        self.enr_velocity_current[:] = self.enr_velocity_new[:]
        self.enr_coordinates_current[:] = self.enr_coordinates_new[:]
        
        # Cohesive model
        self.cohesive_force.increment_values()
        self.damage_variable.increment_values()
        self.discontinuity_opening.increment_values()
        self.dissipated_energy.increment_values()
