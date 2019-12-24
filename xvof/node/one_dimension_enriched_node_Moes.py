#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining one dimension enriched nodes
"""
import numpy as np
from xvof.node.one_dimension_node import OneDimensionNode
from xvof.node.one_dimension_enriched_node import OneDimensionEnrichedNode
from xvof.discontinuity.discontinuity import Discontinuity


class OneDimensionMoesEnrichedNode(OneDimensionEnrichedNode):

    def __init__(self, nbr_of_nodes, initial_positions, initial_velocities, section=1.):
        super(OneDimensionMoesEnrichedNode, self).__init__(nbr_of_nodes, initial_positions,
                                                          initial_velocities, section=section)
        self._v_field = np.copy(self._upundemi)

    @property
    def enrichment_concerned(self):
        """
        :return: boolean mask indicating which nodes are concerned by enrichment
        (fonctions de formes enrichies N* non nulles sur leur support)
        en 1D : correspond aux noeuds enrichis et les premiers voisins.
        """
        # implémentation pour une unique discontinuité,avec 2 noeuds enrichis.
        enrichment_etendu = ~ self.classical
        if np.any(self.enriched):
            index0 = self.enriched.nonzero()[0][0]
            if not self.enriched[index0-1]:
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
        return self._v_field

    def compute_complete_velocity_field(self):
        """
        Calcule le champ de vitesse vrai à partir des champs classiques et enrichis
        """
        self._v_field = np.copy(self._upundemi)
        for disc in Discontinuity.discontinuity_list():
            node_in = np.where(disc.mask_in_nodes)[0][0]
            node_out = np.where(disc.mask_out_nodes)[0][0]
            self._v_field[node_in] -= disc._additional_dof_velocity_new[0, ]  # - u1* à t n+1/2
            self._v_field[node_out] += disc._additional_dof_velocity_new[1, ]  # + u2* à t n+1/2

    def infos(self, index):
        """
        Print information
        """
        OneDimensionNode.infos(self, index)
        message = "==> classical velocity at t-1/2 = {}\n". \
            format(self.umundemi[index])
        message += "==> classical velocity at t+1/2 = {}\n". \
            format(self.upundemi[index])
        message += "==> classical force = {}\n". \
            format(self.force[index])
        message += "------------"
        for disc in [d for d in Discontinuity.discontinuity_list() if index in [d.mask_in_nodes, d.mask_out_nodes]]:
            recal_ind = (disc.mask_out_nodes == index)
            message += "==> enriched velocity at t-1/2 = {}\n". \
                format(disc.additional_dof_velocity_current[recal_ind])
            message += "==> enriched velocity at t+1/2 = {}\n". \
                format(disc.additional_dof_velocity_current[recal_ind])
            message += "==> enriched force = {}\n". \
                format(disc.additional_dof_force[recal_ind])
        print message

    def coupled_enrichment_terms_compute_new_velocity(self, delta_t, inv_matrice_couplage):
        """
        Compute the coupled terms between classical and enriched dof due to non diagonal complete mass matrix
        Takes into account nodes concerned by enrichment and not only the enriched nodes
        :param delta_t: time step
        :param inv_matrice_couplage : inverse de la matrice de masse partie couplage ddl classique / enrichis
        """
        for disc in Discontinuity.discontinuity_list():
            # creation artificielle d'un mask "enrichment_concerned" pour la discontinuity disc
            node_in = np.where(disc.mask_in_nodes)[0][0]
            node_out = np.where(disc.mask_out_nodes)[0][0]
            mask_disc = [node_in - 1, node_in, node_out, node_out + 1]
            # correction de U1* et U2*
            disc._additional_dof_velocity_new += np.dot(
                inv_matrice_couplage.transpose(), self._force[mask_disc]) * delta_t
            # correction de U0 U1 U2 U3 après couplage
            self._upundemi[mask_disc] += np.dot(inv_matrice_couplage, disc.additional_dof_force) * delta_t

    def enriched_nodes_compute_new_coordinates(self, delta_t):
        """
        Compute the new nodes coordinates after enrichment
        :param delta_t: float, time step
        """
        for disc in Discontinuity.discontinuity_list():
            self._xtpdt[disc.mask_in_nodes] -= disc.additional_dof_velocity_new[0] * delta_t
            self._xtpdt[disc.mask_out_nodes] += disc.additional_dof_velocity_new[1] * delta_t

    def enriched_nodes_compute_new_force(self, topology, vecteur_pression_classique, vecteur_pseudo_classique):
        """
        Compute the enriched force on enriched nodes and apply correction for classical force on enriched nodes
        (classical ddl)
        :param topology: Topology1D, give the connectivity of nodes
        :param vecteur_pression_classique: array1D
        :param vecteur_pseudo_classique: array1D
        """
        connectivity = topology.cells_in_contact_with_node[self.enriched].flatten()
        connectivity = np.unique(connectivity)  # ne garder que les n°des cell proches de l'enrichissement sans doublons
        for disc in Discontinuity.discontinuity_list():
            vecteur_pression_enrichie = disc.additional_dof_pressure.new_value
            vecteur_pseudo_enrichie = disc.additional_dof_artificial_viscosity.new_value

            p_enr = vecteur_pression_enrichie + vecteur_pseudo_enrichie  # scalaire
            p_classic = vecteur_pression_classique[connectivity] + vecteur_pseudo_classique[connectivity]
            # vecteur de taille 3
            alpha = 2. * disc.position_in_ruptured_element - 1.

            # Node left
            disc.additional_dof_force[0] = (- p_classic[0] + alpha * p_classic[1] - p_enr) * self.section
            self._force[disc.mask_in_nodes] += alpha * p_enr * self.section

            # Node right
            disc.additional_dof_force[1] = (- alpha * p_classic[1] + p_enr - p_classic[2]) * self.section
            self._force[disc.mask_out_nodes] -= alpha * p_enr * self.section
