# -*- coding: iso-8859-1 -*-
"""
Implementing EnrichElement class
"""
import numpy as np

from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.rupturetreatment.rupturetreatment import RuptureTreatment
from xfv.src.data.enriched_mass_matrix_props import EnrichedMassMatrixProps


class EnrichElement(RuptureTreatment):
    """
    A treatment that enrich one of the ruptured cells
    """
    __never_enriched = True
    __debug = False

    def __init__(self, position_rupture: float, lump_matrix: EnrichedMassMatrixProps):
        super(EnrichElement, self).__init__()
        self.__position_rupture = position_rupture
        self.__lump = lump_matrix

    @property
    def position_rupture(self) -> float:
        """
        Accessor on the relative position of discontinuity inside the enriched cell
        """
        return self.__position_rupture

    @property
    def lump_style(self) -> EnrichedMassMatrixProps:
        """
        Accessor on the mass matrix lumping to be applied
        """
        return self.__lump

    def apply_treatment(self, cells, ruptured_cells, nodes, topology, time):
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
                print("Cells checking rupture criterion are :", np.nonzero(ruptured_cells))
                print("Cells already enriched are :", np.nonzero(cells.enriched))
                print("New cells to be enriched are :", np.nonzero(cells_to_be_enr))

            for cell_tb_enr in np.nonzero(cells_to_be_enr)[0]:
                if not cells.enriched[cell_tb_enr]:
                    print("---------------------------------------------")
                    print("New ruptured cell detected. "
                          "Starting enrichment process for cell {:}".format(cell_tb_enr))
                    print("Beginning enrichment sequence at time {:}".format(time))
                    print("==> Enrichment of cell : ", cell_tb_enr)
                    # Identify nodes to be enriched
                    nodes_to_be_enr = topology.nodes_belonging_to_cell[cell_tb_enr]
                    in_nodes = np.zeros(nodes.number_of_nodes, dtype=bool, order="C")
                    out_nodes = np.zeros(nodes.number_of_nodes, dtype=bool, order="C")
                    in_nodes[nodes_to_be_enr[0]] = True
                    out_nodes[nodes_to_be_enr[1]] = True
                    print("==> Enrichment of nodes : ", nodes_to_be_enr.flatten())
                    print("==> In nodes : ", np.nonzero(in_nodes))
                    print("==> Out nodes : ", np.nonzero(out_nodes))
                    nodes.classical[nodes_to_be_enr] = False
                    # Calcul coordonnées de la disc
                    x_left = nodes.xt[nodes_to_be_enr[0]]
                    x_right = nodes.xt[nodes_to_be_enr[1]]
                    assert x_left < x_right
                    d_coord = x_left + (x_right - x_left) * self.__position_rupture
                    print("==> Discontinuity position : ", x_left, " < ", d_coord, " < ", x_right)
                    # Build the discontinuity
                    disc = Discontinuity(cell_tb_enr, in_nodes, out_nodes,
                                         self.__position_rupture, self.__lump)
                    cells.classical[cell_tb_enr] = False
                    # Initialisation de la partie droite des champs + cell size
                    self.initialize_cracked_cell_size(cells, cell_tb_enr)
                    cells.initialize_additional_cell_dof(disc)
                    nodes.initialize_additional_node_dof(disc)
                else:
                    raise NotImplementedError("""Cell {:} is already enriched.
                    Impossible to enrich it twice""". format(cell_tb_enr))

        ruptured_cells[:] = False

    def initialize_cracked_cell_size(self, cells, cell_tb_enr):
        """
        Compute the size of the each part of the newly enriched cell
        :param cells: cell collection
        :param cell_tb_enr: id if the cell to be enriched
        :return:
        """
        # Affectation left / right part sizes
        cells.right_part_size.new_value = \
            (1. - self.__position_rupture) * cells.size_t_plus_dt[cell_tb_enr]
        cells.left_part_size.new_value = \
            self.__position_rupture * cells.size_t_plus_dt[cell_tb_enr]
        # L'initialisation des tailles gauches et droites courantes n'est
        # pas nécessaire. On initialise simplement avec des tailles fictives de
        # sorte qu'on peut gérer le calcul des forces cohésives à l'itération où la
        # discontinuité est créée. Cette taille ficitive permet simplement
        # d'obtenir une ouverture nulle de la fissure à l'itération de création de
        # la discontinuité.  Elle sera écrasée après ce calcul lors de l'appel
        # de mesh.increment().
        cells.right_part_size.current_value = \
            (1. - self.__position_rupture) * cells.size_t[cell_tb_enr]
        cells.left_part_size.current_value = \
            self.__position_rupture * cells.size_t[cell_tb_enr]