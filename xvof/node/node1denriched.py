#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe dÃ©finissant un noeud enrichi en 1d
"""
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ########### IMPORTATIONS DIVERSES  ####################
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import numpy as np
from xvof.node import Node1d


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ###### DEFINITION DES CLASSES & FONCTIONS  ###############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class Node1dEnriched(Node1d):
    """
    Une classe pour les noeuds enrichis dans le cas 1d

    @todo :
    - Faut il recalculer la masse associÃ©e au noeud en cas d'enrichissement quand
    la discontinuitÃ© n'est pas au milieu de l'Ã©lÃ©ment?
    """
    # pylint: disable-msg=R0902
    # 9 attributs : cela semble raisonnable pour ce cas
    def __init__(self, origin_node, mask, pos_disc):
        '''
        :param origin_node: l'ensemble des noeuds d'origine
        :type origin_node: Node1d
        :param mask: liste des indices des noeuds enrichis
        :type mask: numpy.array
        :param pos_disc: position absolue de la rupture
        :type pos_disc: float
        '''
        nbr_of_nodes = origin_node.number_of_nodes
        poz_init = origin_node.coordt[:]
        vit_init = origin_node.umundemi[:]
        sect = origin_node.section
        Node1d.__init__(self, nbr_of_nodes, poz_init, vit_init, section=sect)
        self._mask = mask
        self._mask_in_nodes = self._xt[mask] - pos_disc < 0  # Noeuds à gauche de la rupture (in)
        self._mask_out_nodes = self._xt[mask] - pos_disc > 0 # Noeuds à droite de la rupture (out)
        #
        self._xtpdt = origin_node.coordtpdt[:]
        self._upundemi = origin_node.upundemi[:]
        #
        self._umundemi_classique = origin_node.umundemi[:]
        self._upundemi_classique = origin_node.upundemi[:]
        self._force_classique = origin_node.force[:]
        # ==> Toutes les variables enrichies sont initialisÃ©es Ã  0
        self._umundemi_enrichi = np.zeros([self.number_of_nodes, self.dimension], dtype=np.float64, order='C')
        self._upundemi_enrichi = np.zeros([self.number_of_nodes, self.dimension], dtype=np.float64, order='C')
        self._force_enrichi = np.zeros([self.number_of_nodes, 1], dtype=np.float64, order='C')
        #
        self._masse = origin_node.masse[:]
        self._invmasse = origin_node.invmasse[:]

    # ------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    # ------------------------------------------------------------
    @property
    def umundemi_classique(self):
        """
        Vitesse classique au demi pas de temps prÃ©cÃ©dent
        """
        return self._umundemi_classique

    @property
    def umundemi_enrichi(self):
        """
        Vitesse enrichie au demi pas de temps prÃ©cÃ©dent
        """
        return self._umundemi_enrichi

    @property
    def upundemi_classique(self):
        """
        Vitesse classique au demi pas de temps suivant
        """
        return self._upundemi_classique

    @property
    def upundemi_enrichi(self):
        """
        Vitesse enrichie au demi pas de temps suivant
        """
        return self._upundemi_enrichi

    @property
    def force_classique(self):
        """
        Force classique
        """
        return self._force_classique

    @property
    def force_enrichi(self):
        """
        Force enrichie
        """
        return self._force_enrichi

    # ------------------------------------------------------------
    # DEFINITIONS DES METHODES
    # ------------------------------------------------------------

    def infos(self, index):
        """
        Affichage des informations
        """
        Node1d.infos(self, index)
        message = "==> vitesse classique Ã  t-1/2 = {}\n".\
            format(self.umundemi_classique[index])
        message += "==> vitesse enrichie Ã  t-1/2 = {}\n".\
            format(self.umundemi_enrichi[index])
        message += "==> vitesse classique Ã  t+1/2 = {}\n".\
            format(self.upundemi_classique[index])
        message += "==> vitesse enrichie Ã  t+1/2 = {}\n".\
            format(self.upundemi_enrichi[index])
        message += "==> force classique = {}\n".\
            format(self.force_classique[index])
        message += "==> force enrichie = {}\n".\
            format(self.force_enrichi[index])
        print message

    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps supÃ©rieur
        """
        self._upundemi_enrichi[self._mask] = \
            self.force_enrichi[self._mask] * self.invmasse[self._mask] * delta_t + self.umundemi_enrichi[self._mask]
        self._upundemi_classique[self._mask] = \
            self.force_classique[self._mask] * self.invmasse[self._mask] * delta_t + self.umundemi_classique[self._mask]
        self._upundemi[self._mask_in_nodes] = \
            self.upundemi_classique[self._mask_in_nodes] - self.upundemi_enrichi[self._mask_in_nodes]
        self._upundemi[self._mask_out_nodes] = \
            self.upundemi_classique[self._mask_out_nodes] + self.upundemi_enrichi[self._mask_out_nodes]

    def calculer_nouvo_force(self, topologie, vecteur_pression_classique, vecteur_pseudo_classique,
                             vecteur_pression_enrichie, vecteur_pseudo_enrichie):
        """
        Calcul des forces agissant sur les noeuds

        :param topologie: topologie du calcul
        :param vecteur_pression_maille: vecteur des pressions de chaque élément + 2 pressions nulles à gauche et à droite
        :param vecteur_pseudo_maille: vecteur des pseudoviscosite de chaque élément

        :type topologie: Topology
        :type vecteur_pression_maille: numpy.array([nbr_of_nodes+2, 1], dtype=np.float64, order='C')
        :type vecteur_pseudo_maille: numpy.array([nbr_of_nodes, 1], dtype=np.int64, order='C')
        """
        # Suppose les éléments voisins triés par position croissante
        connectivity = np.array(topologie._cells_in_contact_with_node[1:self.number_of_nodes - 1])
        # Restriction sur les éléments concernés par l'enrichissement
        connectivity_in = connectivity[self._mask_in_nodes]
        connectivity_out = connectivity[self._mask_out_nodes]
        p_classic = vecteur_pression_classique[connectivity_out] + vecteur_pseudo_classique[connectivity_out]
        p_enr = vecteur_pression_enrichie[connectivity_out] + vecteur_pseudo_enrichie[connectivity_out]
        if self._lev == 1:
            # Noeud Ã  droite de la discontinuitÃ©
            pgauche = elements_voisins[0].pressure.classical_part.new_value + \
                elements_voisins[0].pseudo.classical_part.current_value
            pgauche_enr = \
                elements_voisins[0].pressure.enriched_part.new_value + \
                elements_voisins[0].pseudo.enriched_part.current_value
            pdroite = elements_voisins[1].pressure.new_value + \
                elements_voisins[1].pseudo.current_value
            #
            self._force_classique[:] = (pgauche - pdroite) * self.section
            self._force_enrichi[:] = (pgauche_enr - pdroite) * self.section
        elif self.position_relative == -1:
            # Noeud Ã  gauche de la discontinuitÃ©
            pgauche = elements_voisins[0].pressure.new_value + \
                elements_voisins[0].pseudo.current_value
            pdroite = elements_voisins[1].pressure.classical_part.new_value + \
                elements_voisins[1].pseudo.classical_part.current_value
            pdroite_enr = \
                elements_voisins[1].pressure.enriched_part.new_value + \
                elements_voisins[1].pseudo.enriched_part.current_value
            #
            self._force_classique[:] = (pgauche - pdroite) * self.section
            self._force_enrichi[:] = (-pgauche - pdroite_enr) * self.section
        self._force = None

    def incrementer(self):
        """
        Mise Ã  jour de la vitesse et de la coordonnÃ©e du noeud
        pour passer au pas de temps suivant.
        """
        Node1d.incrementer(self)
        self._umundemi_classique[:] = self.upundemi_classique[:]
        self._umundemi_enrichi[:] = self.upundemi_enrichi[:]

if __name__ == "__main__":
    NODE_INI = Node1d(123, section=1.0e-06)
    MY_NODE = Node1dEnriched(NODE_INI)
    MY_NODE.infos()
