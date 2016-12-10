#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining one dimension enriched nodes
"""
import numpy as np

from xvof.node.one_dimension_node import OneDimensionNode


class OneDimensionEnrichedNode(OneDimensionNode):
    """
    A class for enriched nodes in 1d case.

    @todo :
    - Do we need to compute the mass of the node in case of enrichment when discontinuity is not at the middle of
    the cell?
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
        # self._classical indicate which nodes are classical (not enriched)
        self._classical = np.empty([self.number_of_nodes, ], dtype=np.bool, order='C')
        self._classical[:] = True
        # ==> All enriched field are initialized with a nill value
        self._umundemi_enriched = np.zeros([self.number_of_nodes, self.dimension], dtype=np.float64, order='C')
        self._upundemi_enriched = np.zeros([self.number_of_nodes, self.dimension], dtype=np.float64, order='C')
        self._force_enriched = np.zeros([self.number_of_nodes, 1], dtype=np.float64, order='C')
        # A list of discontinuity objects
        self.discontinuities = []

    @property
    def umundemi_enriched(self):
        """
        Vitesse enrichie au demi pas de temps précédent
        """
        return self._umundemi_enriched

    @property
    def upundemi_enriched(self):
        """
        Vitesse enrichie au demi pas de temps suivant
        """
        return self._upundemi_enriched

    @property
    def force_enriched(self):
        """
        Force enrichie
        """
        return self._force_enriched

    @property
    def classical(self):
        """
        :return: boolean mask indicating which nodes are classical
        """
        return self._classical
    
    @property
    def enriched(self):
        """
        self.enriched indisque les noeus enrichis
        """
        return ~self._classical

    @property
    def velocity_field(self):
        """
        Champ de vitesse vraie
        """
        res = self.upundemi
        for disc in self.discontinuities:
            res[disc.mask_in_nodes] -= self.upundemi_enriched[disc.mask_in_nodes]
            res[disc.mask_out_nodes] += self.upundemi_enriched[disc.mask_out_nodes]
        return res

    def infos(self, index):
        """
        Affichage des informations
        """
        OneDimensionNode.infos(self, index)
        message = "==> classical velocity at t-1/2 = {}\n". \
            format(self.umundemi[index])
        message += "==> enriched velocity at t-1/2 = {}\n". \
            format(self.umundemi_enriched[index])
        message += "==> classical velocity at t+1/2 = {}\n". \
            format(self.upundemi[index])
        message += "==> enriched velocity at t+1/2 = {}\n". \
            format(self.upundemi_enriched[index])
        message += "==> classical force = {}\n". \
            format(self.force[index])
        message += "==> enriched force = {}\n". \
            format(self.force_enriched[index])
        print message

    def enriched_nodes_compute_new_velocity(self, delta_t):
        # Calcul du vecteur vitesse enrichie des noeuds enrichis
        self._upundemi_enriched[self.enriched] = self.force_enriched[self.enriched] * \
                                                 self.invmasse[self.enriched] * delta_t + \
                                                 self.umundemi_enriched[self.enriched]

    def enriched_nodes_compute_new_coordinates(self, delta_t):
        for disc in self.discontinuities:
            self._xtpdt[disc.mask_in_nodes] -= self.upundemi_enriched[disc.mask_in_nodes] * delta_t
            self._xtpdt[disc.mask_out_nodes] += self.upundemi_enriched[disc.mask_out_nodes] * delta_t

    def enriched_nodes_compute_new_force(self, topology, vecteur_pression_classique, vecteur_pression_enrichie,
                                         vecteur_pseudo_classique, vecteur_pseudo_enrichie):
        connectivity = topology.cells_in_contact_with_node[1:-1]
        for disc in self.discontinuities:
            connectivity_in = connectivity[disc.mask_in_nodes[1:-1]].flatten()
            connectivity_out = connectivity[disc.mask_out_nodes[1:-1]].flatten()
            p_classic = vecteur_pression_classique[connectivity_in] + vecteur_pseudo_classique[connectivity_in]
            p_enr = vecteur_pression_enrichie[connectivity_in] + vecteur_pseudo_enrichie[connectivity_in]
            self._force_enriched[disc.mask_in_nodes[1:-1]] = (- p_classic[0] - p_enr[1]) * self.section
            p_classic = vecteur_pression_classique[connectivity_out] + vecteur_pseudo_classique[connectivity_out]
            p_enr = vecteur_pression_enrichie[connectivity_out] + vecteur_pseudo_enrichie[connectivity_out]
            self._force_enriched[disc.mask_out_nodes[1:-1]] = (p_enr[0] - p_classic[1]) * self.section

    def enriched_nodes_increment(self):
        """
        Mise à jour de la vitesse et de la coordonnée du noeud pour passer au pas de temps suivant.
        """
        self._umundemi_enriched[:] = self.upundemi_enriched[:]
