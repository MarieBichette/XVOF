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

    def applyTreatment(self, cells, ruptured_cells, nodes, topology):
        """
        Apply the rupture treatment by enriching one of the cells that is marked as ruptured cells

        :param cells: array of all cells
        :param ruptured_cells: boolean array marking the ruptured cells
        :param nodes: array of all nodes
        :param topology: topology of the problem
        """
        if ruptured_cells.any() and not cells.enriched.any():  # Enrichment is made once for all
            cells_to_be_enr = ruptured_cells
            # We keep only one cell to be enriched
            indices_cells_to_be_enr = np.where(cells_to_be_enr == True)
            cells_to_be_enr[:] = False
            cells_to_be_enr[indices_cells_to_be_enr[0][0]] = True
            #
            nodes_to_be_enr = np.array(topology._nodes_belonging_to_cell)[ruptured_cells]
            print "==> Enrichment of nodes : ", nodes_to_be_enr
            nodes._classiques[nodes_to_be_enr] = False
            nodes.pos_disc = np.mean(nodes.xt[nodes_to_be_enr])
            print "==> Enrichement of cells : ", np.where(cells_to_be_enr == True)
            cells.classical[cells_to_be_enr] = False
            cells.right_size.new_value = (1. - self.__position_rupture) * cells.size_t_plus_dt
            cells.left_size.new_value = self.__position_rupture * cells.size_t_plus_dt
        ruptured_cells[:] = False
