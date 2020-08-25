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
                    nodes_to_be_enr = topology.nodes_belonging_to_cell[cell_tb_enr]
                    nodes_to_be_enr_mask = np.ndarray(nodes.number_of_nodes,
                                                      dtype=np.bool, order="C")
                    nodes_to_be_enr_mask[:] = False
                    nodes_to_be_enr_mask[nodes_to_be_enr.flatten()] = True
                    print("==> Enrichment of nodes : ", nodes_to_be_enr.flatten())
                    nodes.classical[nodes_to_be_enr] = False
                    x_left_node = nodes.xt[nodes_to_be_enr].flatten()[0]
                    x_right_node = nodes.xt[nodes_to_be_enr].flatten()[1]
                    discontinuity_coord = x_left_node + \
                                          (x_right_node - x_left_node) * self.__position_rupture
                    print("==> Discontinuity position : ",
                          x_left_node, " < ", discontinuity_coord, " < ", x_right_node)
                    in_nodes = np.logical_and(
                        nodes_to_be_enr_mask, nodes.xt[:, 0] < discontinuity_coord)
                    out_nodes = np.logical_and(
                        nodes_to_be_enr_mask, nodes.xt[:, 0] > discontinuity_coord)
                    print("==> In nodes : ", np.nonzero(in_nodes))
                    print("==> Out nodes : ", np.nonzero(out_nodes))
                    # Build the discontinuity
                    disc = Discontinuity(in_nodes, out_nodes, self.__position_rupture, self.__lump)
                    disc.find_ruptured_cell_id(topology)
                    cells.classical[cell_tb_enr] = False
                    # Affectation left / right part sizes
                    cells.right_part_size.new_value = \
                        (1. - self.__position_rupture) * cells.size_t_plus_dt[cell_tb_enr]
                    cells.left_part_size.new_value = \
                        self.__position_rupture * cells.size_t_plus_dt[cell_tb_enr]
                    # L'initialisation des tailles gauches et droites courantes n'est
                    # pas n�cessaire. On initialise simplement avec des tailles fictives de
                    # sorte qu'on peut g�rer le calcul des forces coh�sives � l'it�ration o� la
                    # discontinuti� est cr��e. Cette taille ficitive permet simplement
                    # d'obtenir une ouverture nulle de la fissure � l'it�ration de cr�ation de
                    # la discontinuit�.  Elle sera �cras�e apr�s ce calcul lors de l'appel
                    # de mesh.increment().
                    cells.right_part_size.current_value = \
                        (1. - self.__position_rupture) * cells.size_t[cell_tb_enr]
                    cells.left_part_size.current_value = \
                        self.__position_rupture * cells.size_t[cell_tb_enr]
                    # Initialisation de la partie droite des champs pour Hansbo method
                    if not disc.initialisation:
                        cells.initialize_additional_cell_dof(disc)
                        nodes.initialize_additional_node_dof(disc)
                        disc.have_dof_been_initialized()

                    print("---------------------------------------------")
                else:
                    raise NotImplementedError("""La cell {:} est d�j� enrichie.
                    On ne peut pas l'enrichir plusieurs fois""". format(cell_tb_enr))

        ruptured_cells[:] = False
