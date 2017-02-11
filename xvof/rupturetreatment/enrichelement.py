# -*- coding: iso-8859-1 -*-
"""
Implementing EnrichElement class
"""
import numpy as np

from xvof.discontinuity.discontinuity import Discontinuity
from xvof.rupturetreatment.rupturetreatment import RuptureTreatment


class EnrichElement(RuptureTreatment):
    """
    A treatment that enrich one of the ruptured cells
    """
    __never_enriched = True
    __debug = False

    def __init__(self, position_rupture):
        super(EnrichElement, self).__init__()
        self.__position_rupture = position_rupture

    @property
    def position_rupture(self):
        return self.__position_rupture


    def applyTreatment(self, cells, ruptured_cells, nodes, topology):
        """
        Apply the rupture treatment by enriching one of the cells that is marked as ruptured cells

        :param cells: array of all cells
        :param ruptured_cells: boolean array marking the ruptured cells
        :param nodes: array of all nodes
        :param topology: topology of the problem
        """
        if ruptured_cells.any():  # Enrichment is made once for all
            cells_to_be_enr = np.logical_and(ruptured_cells, ~cells.enriched)
            if EnrichElement.__debug:
                print "Beginning enrichment sequence"
                print "Cells checking rupture criterion are :", np.nonzero(ruptured_cells)
                print "Cells already enriched are :", np.nonzero(cells.enriched)
                print "New cells to be enriched are :", np.nonzero(cells_to_be_enr)
            for cell_tb_enr in np.nonzero(cells_to_be_enr)[0]:
                if Discontinuity.discontinuity_number() < 1:
                    print "==> Enrichment of cell : ", cell_tb_enr
                    nodes_to_be_enr = topology.nodes_belonging_to_cell[cell_tb_enr]
                    nodes_to_be_enr_mask = np.ndarray(nodes.number_of_nodes, dtype=np.bool, order="C")
                    nodes_to_be_enr_mask[:] = False
                    nodes_to_be_enr_mask[nodes_to_be_enr.flatten()] = True
                    print "==> Enrichment of nodes : ", nodes_to_be_enr.flatten()
                    nodes.classical[nodes_to_be_enr] = False
                    x_left_node = nodes.xt[nodes_to_be_enr].flatten()[0]
                    x_right_node = nodes.xt[nodes_to_be_enr].flatten()[1]
                    discontinuity_coord = x_left_node + (x_right_node - x_left_node) * self.__position_rupture
                    print "==> Discontinuity position : ", x_left_node, " < ", discontinuity_coord, " < ", x_right_node
                    in_nodes = np.logical_and(nodes_to_be_enr_mask, nodes.xt[:, 0] < discontinuity_coord)
                    out_nodes = np.logical_and(nodes_to_be_enr_mask, nodes.xt[:, 0] > discontinuity_coord)
                    print "==> In nodes : ", np.nonzero(in_nodes)
                    print "==> Out nodes : ", np.nonzero(out_nodes)
                    Discontinuity.discontinuity_list().append(Discontinuity(in_nodes, out_nodes))
                    cells.classical[cell_tb_enr] = False
                    cells.right_size.new_value = (1. - self.__position_rupture) * cells.size_t_plus_dt
                    cells.left_size.new_value = self.__position_rupture * cells.size_t_plus_dt
        ruptured_cells[:] = False
