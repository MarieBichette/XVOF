# -*- coding: iso-8859-1 -*-
"""
Implementing the MassMatrix class
"""
import numpy as np


def compute_wilkins_mass_matrix(topology, cell_mass_vector, node_number_by_cell_vector):
    """
    Compute nodal mass by averaging the mass of neighbouring cells (Wilkins method)

    :param topology: topology of the simulation
    :param cell_mass_vector: cells mass vector
    :param node_number_by_cell_vector: number of nodes per cell (vector)

    :type topology: Topology
    :type cell_mass_vector: numpy.array([nbr_of_nodes, 1], dtype=np.float64, order='C')
    :type node_number_by_cell_vector: numpy.array([nbr_of_nodes, 1], dtype=np.int64, order='C')
    """
    nbr_nodes = len(topology.cells_in_contact_with_node)
    mass_matrix = np.zeros([nbr_nodes, 1], dtype=np.float64, order='C')
    for ind_node in xrange(nbr_nodes):
        neighboring_cells = topology.getCellsInContactWithNode(ind_node)
        # Index -1 => unexistant cell
        neighboring_cells = neighboring_cells[neighboring_cells != -1]
        mass_matrix[ind_node] = np.sum(
            cell_mass_vector[neighboring_cells] / node_number_by_cell_vector[neighboring_cells])
    return mass_matrix
