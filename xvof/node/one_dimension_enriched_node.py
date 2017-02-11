#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining one dimension enriched nodes
"""
import numpy as np

from xvof.node.one_dimension_node import OneDimensionNode
from xvof.discontinuity.discontinuity import discontinuity_list

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
        :return: boolean mask indicating which nodes are enriched
        """
        return ~ self.classical

    @property
    def enrichment_concerned(self):
        """
        :return: boolean mask indicating which nodes are concerned by enrichment (fonctions de formes enrichies N* non nulles sur leur support)
        en 1D : correspond aux noeuds enrichis et les premiers voisins.
        """
        #implémentation pour une unique discontinuité,avec 2 noeuds enrichis.
        enrichment_etendu = ~ self.classical
        if np.any(self.enriched):
            index0 = self.enriched.nonzero()[0][0]
            if self.enriched[index0-1] == False :
                enrichment_etendu[index0-1] = True
                enrichment_etendu[index0 + 2] = True
        return enrichment_etendu

    @property
    def enrichment_not_concerned(self):
        """
        :return: boolean mask indicating which nodes are not concerned by enrichment
        """
        return ~ self.enrichment_concerned

    @property
    def velocity_field(self):
        """
        Champ de vitesse vraie
        """
        return self.compute_complete_velocity_field()

    def compute_complete_velocity_field(self):
        """
        Calcule le champ de vitesse vrai à partir des champs classiques et enrichis
        """
        res = np.copy(self._upundemi)
        for disc in discontinuity_list:
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

    def enriched_nodes_compute_new_velocity(self, delta_t, mask, inv_matrice_masse):
        """
        Compute the new velocity enriched degree of freedom
        :param delta_t: float, time step
        :param mask: array of boolean indicating the location of enriched nodes
        :param matrice_masse: mass matrix to be inverted
        :return:
        """
        # utilise un mask et une matrice masse particulière pour garder la forme de compute_new_velocity classique.
        # permettra de gérer plusieurs discontinuités

        # Ancien : Calcul du vecteur vitesse enrichie des noeuds enrichis
        #self._upundemi_enriched[self.enriched] = \
        # self.force_enriched[self.enriched] * self.invmasse[self.enriched] * delta_t + self.umundemi_enriched[self.enriched]

        self._upundemi_enriched[mask] = \
            np.dot(inv_matrice_masse , self.force_enriched[mask]) * delta_t + self.umundemi_enriched[mask]

    def coupled_enrichment_terms_compute_new_velocity(self, delta_t, inv_matrice_couplage):
        """
        Compute the coupled terms between classical and enriched dof due to non diagonal complete mass matrix
        Takes into account nodes concerned by enrichment and not only the enriched nodes
        :param delta_t: time step
            MARCHE POUR UNE SEULE DISCONTINUITE
        """
        # calcul de U1* et U2*
        self._upundemi_enriched[self.enriched] += np.dot(inv_matrice_couplage.transpose(),
                                                         self.force[self.enrichment_concerned]) * delta_t
        # calcul de U0 U1 U2 U3 après couplage
        self._upundemi[self.enrichment_concerned] += np.dot(inv_matrice_couplage,
                                                            self.force_enriched[self.enriched]) * delta_t

    def enriched_nodes_compute_new_coordinates(self, delta_t):
        """
        Compute the new nodes coordinates after enrichment
        :param delta_t: float, time step
        """
        for disc in discontinuity_list:
            self._xtpdt[disc.mask_in_nodes] -= self.upundemi_enriched[disc.mask_in_nodes] * delta_t
            self._xtpdt[disc.mask_out_nodes] += self.upundemi_enriched[disc.mask_out_nodes] * delta_t

    def enriched_nodes_compute_new_force(self, topology, vecteur_pression_classique, vecteur_pression_enrichie,
                                         vecteur_pseudo_classique, vecteur_pseudo_enrichie):
        """
        Compute the enriched force on enriched nodes
        :param topology: Topology1D, give the connectivity of nodes
        :param vecteur_pression_classique: array1D
        :param vecteur_pression_enrichie: array1D
        :param vecteur_pseudo_classique: array1D
        :param vecteur_pseudo_enrichie: array1D
        """
        connectivity = topology.cells_in_contact_with_node[1:-1]
        for disc in discontinuity_list:
            connectivity_in = connectivity[disc.mask_in_nodes[1:-1]].flatten()
            connectivity_out = connectivity[disc.mask_out_nodes[1:-1]].flatten()
            p_classic = vecteur_pression_classique[connectivity_in] + vecteur_pseudo_classique[connectivity_in]
            p_enr = vecteur_pression_enrichie[connectivity_in] + vecteur_pseudo_enrichie[connectivity_in]
            self._force_enriched[disc.mask_in_nodes] = (- p_classic[0] - p_enr[1]) * self.section
            p_classic = vecteur_pression_classique[connectivity_out] + vecteur_pseudo_classique[connectivity_out]
            p_enr = vecteur_pression_enrichie[connectivity_out] + vecteur_pseudo_enrichie[connectivity_out]
            self._force_enriched[disc.mask_out_nodes] = (p_enr[0] - p_classic[1]) * self.section

    def enriched_nodes_increment(self):
        """
        Mise à jour de la vitesse et de la coordonnée du noeud pour passer au pas de temps suivant.
        """
        # self._umundemi_enriched[:] = self.upundemi_enriched[:]
        self._umundemi_enriched = np.copy(self.upundemi_enriched)