# -*- coding: utf-8 -*-
"""
A module implementing the Discontinuity class
"""
import numpy as np
from xfv.src.mesh.topology1d import Topology1D
from xfv.src.fields.field import Field
from xfv.src.data.enriched_mass_matrix_props import EnrichedMassMatrixProps


class Discontinuity:
    """
    A class describing a discontinuity 1D
    """

    # A list of discontinuities
    __discontinuity_list = []
    additional_dof_velocity_current = np.zeros([], dtype=float)
    additional_dof_velocity_new = np.zeros([], dtype=float)
    additional_dof_force = np.zeros([], dtype=float)

    def __init__(self, mask_in_nodes: np.array, mask_out_nodes: np.array,
                 discontinuity_position_in_ruptured_element: float,
                 enriched_mass_matrix_props: EnrichedMassMatrixProps):
        """
        Initializing a single discontinuity after enrichment.
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
            Discontinuity.additional_dof_velocity_current = np.copy(init)
            Discontinuity.additional_dof_velocity_new = np.copy(init)
            Discontinuity.additional_dof_force = np.copy(init)
        else:
            Discontinuity.additional_dof_velocity_current = np.append(Discontinuity.additional_dof_velocity_current, init, axis=0)
            Discontinuity.additional_dof_velocity_new = np.append(Discontinuity.additional_dof_velocity_new, init, axis=0)
            Discontinuity.additional_dof_force = np.append(Discontinuity.additional_dof_force, init, axis=0)
        for ind, disc in enumerate(Discontinuity.discontinuity_list()):
            disc.additional_dof_velocity_current = Discontinuity.additional_dof_velocity_current[ind]
            disc.additional_dof_velocity_new = Discontinuity.additional_dof_velocity_new[ind]
            disc.additional_dof_force = Discontinuity.additional_dof_force[ind]
        self.__mask_in_nodes = mask_in_nodes
        self.__mask_out_nodes = mask_out_nodes

        # Discontinuity cell information
        self._discontinuity_position = discontinuity_position_in_ruptured_element
        self.__mask_ruptured_cell = np.zeros(len(self.mask_in_nodes)-1, dtype=bool)
        self._ruptured_cell_id = 0

        # Indicators of discontinuity state
        self.__mass_matrix_updated = False
        self.__initialisation = False

        # Additional dof representing either the enriched Heaviside value or
        # the field value in the right part of enriched element.
        self._additional_dof_coordinates_current = np.zeros([2, 1])
        self._additional_dof_coordinates_new = np.zeros([2, 1])
        # Damage indicators with cohesive zone model
        # (Always created but null if no damage data in the XDATA.json file...)
        self.cohesive_force = Field(1, current_value=0., new_value=0.)
        self.discontinuity_opening = Field(1, current_value=0., new_value=0.)
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
            if disc._ruptured_cell_id == cell_id:
                return disc
            try_index += 1
        return None

    @property
    def position_in_ruptured_element(self):
        """
        Accessor on the relative position of the discontinuity in ruptured element
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
        Accessor on the boolean that indicates if the mass matrix has been computed
        :return: the boolean that indicates if the mass matrix has been computed
        """
        return self.__mass_matrix_updated

    @property
    def initialisation(self):
        """
        Boolean to follow if the right part values have been correctly initialized for
        Hansbo enrichment. Indeed, those variables have physical meaning and are not zero
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

    def find_ruptured_cell_id(self, topology: Topology1D):
        """
        Compute the cell id of ruptured element associated with this discontinuity
        self.ruptured_cell_id is an integer
        self.mask_ruptured_cell is an array of boolean = True for ruptured cell
        """
        self._ruptured_cell_id = topology.cells_in_contact_with_node[self.mask_in_nodes][0][1]
        self.__mask_ruptured_cell[self._ruptured_cell_id] = True

    def compute_discontinuity_new_opening(self, node_position: np.array):
        """
        Compute the discontinuity opening
        :param node_position: coordinates of the nodes
        """
        epsilon = self._discontinuity_position
        coord_g = node_position[self.mask_in_nodes]
        coord_d = node_position[self.mask_out_nodes]
        enr_coord_g = self._additional_dof_coordinates_new[0]  # x2-
        enr_coord_d = self._additional_dof_coordinates_new[1]  # x1+
        xg_new = (1 - epsilon) * coord_g + epsilon * enr_coord_g
        xd_new = (1 - epsilon) * enr_coord_d + epsilon * coord_d
        self.discontinuity_opening.new_value = (xd_new - xg_new)[0][0]

    @property
    def discontinuity_position(self) -> float:
        """
        Accessor on the discontinuity position
        """
        return self._discontinuity_position

    @property
    def additional_dof_coordinates_current(self):
        """
        Accessor on the additional nodes coordinates at time t
        """
        return self._additional_dof_coordinates_current

    @property
    def additional_dof_coordinates_new(self):
        """
        Accessor on the additional nodes coordinates at time t+dt
        """
        return self._additional_dof_coordinates_new

    def reinitialize_kinematics_after_contact(self):
        """
        Set the new velocity to the old one to cancel the increment that has lead to contact
        """
        Discontinuity.additional_dof_velocity_new[self.__label - 1] = np.copy(Discontinuity.additional_dof_velocity_current[self.__label - 1])
        self._additional_dof_coordinates_new = np.copy(self._additional_dof_coordinates_current)

    def additional_dof_increment(self):
        """
        Increment the variables of discontinuity
        """
        # Kinematics
        self.additional_dof_velocity_current[:] = self.additional_dof_velocity_new[:]
        self._additional_dof_coordinates_current[:] = self.additional_dof_coordinates_new[:]
        # Cohesive model
        self.cohesive_force.increment_values()
        self.damage_variable.increment_values()
        self.discontinuity_opening.increment_values()
