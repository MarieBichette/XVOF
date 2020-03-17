#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining one dimension enriched nodes
"""
import numpy as np

from xvof.src.node.one_dimension_node import OneDimensionNode
from xvof.src.discontinuity.discontinuity import Discontinuity
from xvof.src.mass_matrix.mass_matrix_utilities import multiplicationMasse


class OneDimensionEnrichedNode(OneDimensionNode):
    """
    A class for enriched nodes in 1d case.
    """
    # pylint: disable-msg=R0902
    # 9 attributes : seams of ok here
    def __init__(self, nbr_of_nodes, initial_positions, initial_velocities, section=1.):
        """
        :param nbr_of_nodes: number of nodes
        :type nbr_of_nodes: int
        :param initial_positions: nodes initial positions
        :type initial_positions: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        :param initial_velocities: nodes initial velocities
        :type initial_velocities: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        :param section: section area associated with the node
        :type section: float
        """
        super(OneDimensionEnrichedNode, self).__init__(nbr_of_nodes, initial_positions, initial_velocities,
                                                       section=section)
        self._v_field = np.zeros([self.number_of_nodes])

    @property
    def classical(self):
        """
        commun aux 2
        :return: boolean mask indicating which nodes are classical
        """
        return self._classical
    
    @property
    def enriched(self):
        """
        commun aux 2
        :return: boolean mask indicating which nodes are enriched
        """
        return ~ self.classical

    def compute_additional_dof_new_velocity(self, delta_t, inv_matrice_masse):
        """
        Compute the new velocity enriched degree of freedom
        :param delta_t: float, time step
        :param inv_matrice_masse: inverse of the mass matrix
        :return:
        """
        for disc in Discontinuity.discontinuity_list():
            disc._additional_dof_velocity_new = disc.additional_dof_velocity_current + delta_t * \
                                                     multiplicationMasse(inv_matrice_masse, disc.additional_dof_force)


