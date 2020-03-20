#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining one dimension enriched nodes
"""
import numpy as np
from xvof.src.node.one_dimension_enriched_node import OneDimensionEnrichedNode
from xvof.src.node.one_dimension_node import OneDimensionNode
from xvof.src.discontinuity.discontinuity import Discontinuity
from xvof.src.data.data_container import DataContainer


class OneDimensionHansboEnrichedNode(OneDimensionEnrichedNode):

    def __init__(self, nbr_of_nodes, initial_positions, initial_velocities, section=1.):
        super(OneDimensionHansboEnrichedNode, self).__init__(nbr_of_nodes, initial_positions,
                                                            initial_velocities, section=section)
        self._v_field = np.copy(self._upundemi)

    @property
    def enrichment_concerned(self):
        """
        :return: boolean mask indicating which nodes are concerned by enrichment
        (fonctions de formes enrichies N* non nulles sur leur support)
        en 1D : correspond aux noeuds enrichis et les premiers voisins.
        """
        return ~ self.classical

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
        Calcul du champ de vitesse vraie
        """
        self._v_field = np.copy(self._upundemi)

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
            message += "==> additional degree of freedom for velocity at t-1/2 = {}\n". \
                format(disc.additional_dof_velocity_current[recal_ind])
            message += "==> additional degree of freedom for velocity at t+1/2 = {}\n". \
                format(disc.additional_dof_velocity_current[recal_ind])
            message += "==> additional degree of freedom for force = {}\n". \
                format(disc.additional_dof_force.current_value[recal_ind])  # (pas encore calculé les new forces)
        print message

    def initialize_additional_node_dof(self, disc):
        """
        Initialialise les ddl enrichis aux noeuds
        :param disc:
        :return:
        """
        disc.additional_dof_velocity_current[:] = np.copy(self.umundemi[disc.mask_disc_nodes])
        disc.additional_dof_velocity_new[:] = np.copy(self.upundemi[disc.mask_disc_nodes])

    def enriched_nodes_compute_new_coordinates(self, delta_t):
        """
        Compute the new nodes coordinates after enrichment
        :param delta_t: float, time step
        """
        # Rien ne change dans le calcul de la position des noeuds enrichis.
        # Rappel : les coordonnées enrichies n'existent pas
        OneDimensionNode.compute_new_coodinates(self, delta_t)

    def coupled_enrichment_terms_compute_new_velocity(self, delta_t, inv_matrice_couplage):
        """
        Compute the coupled terms between classical and enriched dof due to non diagonal complete mass matrix
        Takes into account nodes concerned by enrichment and not only the enriched nodes
        :param delta_t: time step
        :param inv_matrice_couplage : inverse de la matrice de masse partie couplage ddl classiq /enrichis
        """
        for disc in Discontinuity.discontinuity_list():
            node_in = np.where(disc.mask_in_nodes)[0][0]
            node_out = np.where(disc.mask_out_nodes)[0][0]
            mask_disc = [node_in, node_out]
            disc._additional_dof_velocity_new += np.dot(inv_matrice_couplage.transpose(),
                                                                self._force[mask_disc]) * delta_t
            self._upundemi[mask_disc] += np.dot(inv_matrice_couplage, disc.additional_dof_force) * delta_t

    def compute_enriched_nodes_new_force(self, topology, contrainte):
        """
        Compute the enriched force on enriched nodes and apply correction for classical force on enriched nodes
        (classical ddl)
        :param topology: Topology1D, give the connectivity of nodes
        :param contrainte : vecteur contrainte xx, array de taille (nb_cell, 1)
        """
        for disc in Discontinuity.discontinuity_list():
            connectivity = topology.getCellsInContactWithNode(topology.getNodesBelongingToCell(disc.ruptured_cell_id))
            assert connectivity.size == 4
            connectivity = connectivity.reshape(2, 2)
            epsilon = disc.position_in_ruptured_element
            node_left_connectivity = connectivity[0]
            node_right_connectivity = connectivity[1]

            assert node_left_connectivity[1] == node_right_connectivity[0]
            cell0 = node_left_connectivity[0]
            cell1 = node_left_connectivity[1]
            cell2 = node_right_connectivity[1]

            sigma0, sigma2 = np.array([0.]), np.array([0.])
            if cell0 != -1:
                sigma0 = contrainte[cell0]
            if cell2 != -1:
                sigma2 = contrainte[cell2]

            sigma1G = contrainte[cell1]
            sigma1D = disc.additional_dof_stress[:,0]

            F1g = sigma1G * epsilon - sigma0
            F2g = - sigma1G * epsilon
            F2d = sigma2 - sigma1D * (1 - epsilon)
            F1d = sigma1D * (1 - epsilon)

            disc.additional_dof_force[0] = F2g * self.section
            disc.additional_dof_force[1] = F1d * self.section
            self._force[disc.mask_in_nodes] = F1g * self.section  # écrase la valeur classique calculée juste avant
            self._force[disc.mask_out_nodes] = F2d * self.section  # écrase la valeur classique calculée juste avant

    def compute_discontinuity_opening(self):
        """
        Compute the opening of discontinuities
        :return:
        """
        for disc in Discontinuity.discontinuity_list():
            # On calcule la nouvelle ouverture de l'écaille
            xd_new = self.xtpdt[disc.mask_out_nodes] - disc.right_part_size.new_value
            xg_new = self.xtpdt[disc.mask_in_nodes] + disc.left_part_size.new_value
            disc.discontinuity_opening.new_value = (xd_new - xg_new)[0][0]

            if (disc.discontinuity_opening.new_value < 0 and
                DataContainer().material_target.damage_model.cohesive_model is None):
                print "Problème avec la discontinuité {:} : ouverture négative".format(disc.label)

    def compute_enriched_nodes_cohesive_forces(self):
        """
        Compute the cohesive forces for the enriched nodes
        """
        self.compute_discontinuity_opening()  # normalement on en a pas besoin mais on le refait au cas où

        for disc in Discontinuity.discontinuity_list():
            # On calcule la nouvelle ouverture de l'écaille
            # xd_new = self.xtpdt[disc.mask_out_nodes] - disc.right_part_size.new_value
            # xg_new = self.xtpdt[disc.mask_in_nodes] + disc.left_part_size.new_value
            # disc.discontinuity_opening.new_value = (xd_new - xg_new)[0][0]

            cohesive_stress = DataContainer().material_target.damage_model.cohesive_model.compute_cohesive_stress(disc)
            disc.cohesive_force.new_value = cohesive_stress

            F_czm = self.section * cohesive_stress

            # On applique la force du ressort sur les forces nodales
            self._force[disc.mask_in_nodes] += (1. - disc.position_in_ruptured_element) * F_czm
            self._force[disc.mask_out_nodes] += - disc.position_in_ruptured_element * F_czm
            disc.additional_dof_force[1] += - (1. - disc.position_in_ruptured_element) * F_czm
            disc.additional_dof_force[0] += disc.position_in_ruptured_element * F_czm



