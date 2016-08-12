# -*- coding: iso-8859-1 -*-
"""
Implementing EnrichElement class
"""
import numpy as np

from xvof.rupturetreatment.rupturetreatment import RuptureTreatment


class EnrichElement(RuptureTreatment):
    """
    A treatment that enrich one of the ruptured cells
    """
    __never_enriched = True

    def __init__(self, position_rupture):
        super(EnrichElement, self).__init__()
        self.__position_rupture = position_rupture

    def applyTreatment(self, cells, ruptured_cells, nodes, topology, cells_coordinates):
        """
        Apply the rupture treatment by enriching one of the cells that is marked as ruptured cells

        :param cells: array of all cells
        :param ruptured_cells: boolean array marking the ruptured cells
        :param nodes: array of all nodes
        :param topology: topology of the problem
        :param cells_coordinates: coordinates of the cells
        """
        cells_to_be_enr = ruptured_cells
        # We keep only one cell to be enriched
        indices_cells_to_be_enr = np.where(cells_to_be_enr == True)
        cells_to_be_enr[:] = False
        cells_to_be_enr[indices_cells_to_be_enr[0][0]] = True
        #
        nodes_to_be_enr = np.array(topology._nodes_belonging_to_cell)[ruptured_cells]
        print "==> Enrichment of nodes : ", nodes_to_be_enr
        nodes._classiques[nodes_to_be_enr] = False
        for pos in cells_coordinates[cells_to_be_enr]:
            nodes.pos_disc = pos[0]
        print "==> Enrichement of cells : ", np.where(cells_to_be_enr == True)
        cells._classical[cells_to_be_enr] = False
        cells._enriched[cells_to_be_enr] = True
        cells.right_size.new_value = self.__position_rupture * cells.size_t_plus_dt
        cells.left_size.new_value = self.__position_rupture * cells.size_t_plus_dt
