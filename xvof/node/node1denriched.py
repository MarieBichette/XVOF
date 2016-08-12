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
    def __init__(self, nbr_of_nodes, poz_init, vit_init, section=1.):
        '''
        :param nbr_of_nodes: nombre de noeuds
        :type nbr_of_nodes: int
        :param poz_init: vecteur position initiale des noeuds
        :type poz_init: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        :param vit_init: vecteur vitesse initiale des noeuds
        :type vit_init: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')      
        '''
        super(Node1dEnriched, self).__init__(nbr_of_nodes, poz_init, vit_init, section=section)
        # self._classiques indique quels noeuds sont classiques (non enrichis)
        self._classiques = np.empty([self.number_of_nodes,], dtype=np.bool, order='C')
        self._classiques[:] = True
        # ==> Toutes les variables enrichies sont initialisÃ©es Ã  0
        self._umundemi_enrichi = np.zeros([self.number_of_nodes, self.dimension], dtype=np.float64, order='C')
        self._upundemi_enrichi = np.zeros([self.number_of_nodes, self.dimension], dtype=np.float64, order='C')
        self._force_enrichi = np.zeros([self.number_of_nodes, 1], dtype=np.float64, order='C')
        # 
        # Un dictionnaire associant à chaque position de discontinuité, les masques des
        # noeuds qui sont in et out
        self._pos_disc = {}

    # ------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    # ------------------------------------------------------------
    @property
    def umundemi_enrichi(self):
        """
        Vitesse enrichie au demi pas de temps précédent
        """
        return self._umundemi_enrichi

    @property
    def upundemi_enrichi(self):
        """
        Vitesse enrichie au demi pas de temps suivant
        """
        return self._upundemi_enrichi

    @property
    def force_enrichi(self):
        """
        Force enrichie
        """
        return self._force_enrichi
    
    @property
    def _enrichis(self):
        '''
        self._enrichis indisque les noeus enrichis
        '''
        return ~self._classiques

    @property
    def velocity_field(self):
        '''
        Champ de vitesse vraie
        '''
        res = self.upundemi
        if self._pos_disc != {}:
            for pos_disc in self.pos_disc.keys():
            # Prise en compte des champs enrichis pour le calcul des nouvelles coordonnées
            # des noeuds enrichis
                mask_in_nodes = self.pos_disc[pos_disc]["inside"]
                mask_out_nodes = self.pos_disc[pos_disc]["outside"]
                res[mask_in_nodes] -= self.upundemi_enrichi[mask_in_nodes]
                res[mask_out_nodes] += self.upundemi_enrichi[mask_out_nodes]
        return res
    # ------------------------------------------------------------
    # DEFINITIONS DES METHODES
    # ------------------------------------------------------------
    @property
    def pos_disc(self):
        return self._pos_disc
    
    @pos_disc.setter
    def pos_disc(self, pos):
        if not self._pos_disc.has_key(pos):
            self._pos_disc[pos] = {}
            mask_in = np.logical_and(self._enrichis, self._xt[:, 0] - pos < 0)
            mask_out = np.logical_and(self._enrichis, self._xt[:, 0] - pos > 0)
            self._pos_disc[pos]["inside"] = mask_in
            self._pos_disc[pos]["outside"] = mask_out

    def infos(self, index):
        """
        Affichage des informations
        """
        Node1d.infos(self, index)
        message = "==> vitesse classique Ã  t-1/2 = {}\n".\
            format(self.umundemi[index])
        message += "==> vitesse enrichie Ã  t-1/2 = {}\n".\
            format(self.umundemi_enrichi[index])
        message += "==> vitesse classique Ã  t+1/2 = {}\n".\
            format(self.upundemi[index])
        message += "==> vitesse enrichie Ã  t+1/2 = {}\n".\
            format(self.upundemi_enrichi[index])
        message += "==> force classique = {}\n".\
            format(self.force_classique[index])
        message += "==> force enrichie = {}\n".\
            format(self.force_enrichi[index])
        print message

    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps supérieur
        """
        # Le vecteur de vitesse des noeuds classique et le vecteur classique des
        # noeuds enrichis sont calculés de la même façon
        super(Node1dEnriched, self).calculer_nouvo_vitesse(delta_t)
        # Calcul du vecteur vitesse enrichie des noeuds enrichis
        self._upundemi_enrichi[self._enrichis] = \
            self.force_enrichi[self._enrichis] * self.invmasse[self._enrichis] * delta_t +\
            self.umundemi_enrichi[self._enrichis]

    def calculer_nouvo_coord(self, delta_t):
        """
        Calcul de la coordonnée au temps t+dt

        :param delta_t: pas de temps
        :type delta_t: float
        """
        # Calcul des nouvelles coordonnées pour les noeuds non enrichis et calcul
        # des nouvelles coordoonées classiques des noeuds enrichis
        super(Node1dEnriched, self).calculer_nouvo_coord(delta_t)
        if self._pos_disc != {}:
            for pos_disc in self.pos_disc.keys():
            # Prise en compte des champs enrichis pour le calcul des nouvelles coordonnées
            # des noeuds enrichis
                mask_in_nodes = self.pos_disc[pos_disc]["inside"]
                mask_out_nodes = self.pos_disc[pos_disc]["outside"]
                self._xtpdt[mask_in_nodes] -= self.upundemi_enrichi[mask_in_nodes] * delta_t
                self._xtpdt[mask_out_nodes] += self.upundemi_enrichi[mask_out_nodes] * delta_t

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
        # Calcul de la force pour les noeuds classiques et de la force classique pour les noeuds
        # enrichis
        super(Node1dEnriched, self).calculer_nouvo_force(topologie, vecteur_pression_classique, vecteur_pseudo_classique)
        if self._pos_disc != {}:
            # Boucle sur les discontinuités
            for pos_disc in self.pos_disc.keys():
                # Suppose les éléments voisins triés par position croissante
                connectivity = np.array(topologie._cells_in_contact_with_node[1:self.number_of_nodes - 1])
                mask_in_nodes = self.pos_disc[pos_disc]["inside"]
                mask_out_nodes = self.pos_disc[pos_disc]["outside"]
                # Restriction sur les éléments concernés par l'enrichissement
                connectivity_out = connectivity[mask_out_nodes][0] 
                p_classic = vecteur_pression_classique[connectivity_out] + vecteur_pseudo_classique[connectivity_out]
                p_enr = vecteur_pression_enrichie[connectivity_out] + vecteur_pseudo_enrichie[connectivity_out]
                self._force_enrichi[mask_out_nodes] = (p_enr[0] - p_classic[1]) * self.section
                connectivity_in = connectivity[mask_in_nodes][0]
                p_classic = vecteur_pression_classique[connectivity_in] + vecteur_pseudo_classique[connectivity_in]
                p_enr = vecteur_pression_enrichie[connectivity_in] + vecteur_pseudo_enrichie[connectivity_in]
                self._force_enrichi[mask_in_nodes] = (- p_classic[0] - p_enr[1]) * self.section

    def incrementer(self):
        """
        Mise Ã  jour de la vitesse et de la coordonnÃ©e du noeud
        pour passer au pas de temps suivant.
        """
        Node1d.incrementer(self)
        self._umundemi_enrichi[:] = self.upundemi_enrichi[:]
