#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un noeud enrichi en 1d
"""
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from xvof.node import Node1d
import numpy as np

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


class Node1dUpgraded(Node1d):
    """
    Une classe pour les noeuds enrichis dans le cas 1d
    """
    # pylint: disable-msg=R0902
    # 9 attributs : cela semble raisonnable pour ce cas
    def __init__(self, origin_node):
        Node1d.__init__(self, origin_node.index, poz_init=origin_node.coordt,
                      vit_init=origin_node.umundemi, section=origin_node.section)

        self._upundemi = origin_node.upundemi[:]
        self._force = np.zeros(1, dtype=float)
        #
        self._umundemi_classique = origin_node.umundemi[:]
        self._upundemi_classique = origin_node.upundemi[:]
        self._force_classique = origin_node.force[:]
        #==> Toutes les variables enrichies sont initialisées à 0
        self._umundemi_enrichi = np.zeros(1, dtype=float)
        self._upundemi_enrichi = np.zeros(1, dtype=float)
        self._force_enrichi = np.zeros(1, dtype=float)
        #
        self.__position_relative = None

    #------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    #------------------------------------------------------------

    @property
    def position_relative(self):
        """
        Position du noeud par rapport à la discontinuité
        """
        return self.__position_relative

    @position_relative.setter
    def position_relative(self, pos):
        """
        Setter de la position relative
        """
        if (pos not in (-1, 1)):
            message = "La position relative du noeud ne peut être que :\n"
            message += " -1 si il est à gauche de la discontinuité\n"
            message += " +1 si il est à droite"
            raise SystemExit(message)
        self.__position_relative = pos

    @property
    def umundemi_classique(self):
        """
        Vitesse classique au demi pas de temps précédent
        """
        return self._umundemi_classique

    @property
    def umundemi_enrichi(self):
        """
        Vitesse enrichie au demi pas de temps précédent
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

    #------------------------------------------------------------
    # DEFINITIONS DES METHODES
    #------------------------------------------------------------

    def infos(self):
        """
        Affichage des informations
        """
        Node1d.infos(self)
        message = "==> vitesse classique à t-1/2 = {}\n".\
        format(self.umundemi_classique)
        message += "==> vitesse enrichie à t-1/2 = {}\n".\
        format(self.umundemi_enrichi)
        message += "==> vitesse classique à t+1/2 = {}\n".\
        format(self.upundemi_classique)
        message += "==> vitesse enrichie à t+1/2 = {}\n".\
        format(self.upundemi_enrichi)
        message += "==> force classique à t-1/2 = {}\n".\
        format(self.force_classique)
        message += "==> force enrichie à t-1/2 = {}\n".\
        format(self.force_enrichi)
        if(self.position_relative is None):
            message += "==> position relative  = None"
        else:
            message += "==> position relative  = {:2d}".\
        format(self.position_relative)
        print message

    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps supérieur

        TEST UNITAIRE
        >>> import numpy as np
        >>> class element:
        ...     pass
        ...
        >>> elem_gauche = element()
        >>> elem_droite = element()
        >>> elem_gauche.coord = np.array([-0.5])
        >>> elem_droite.coord = np.array([0.5])
        >>> elem_gauche.pressure = 2.5e+09
        >>> elem_droite.pressure = 1.0e+09
        >>> elem_gauche.masse = 3./4.
        >>> elem_droite.masse = 1./4.
        >>> node_ini = Node1d(123, section=1.0e-06)
        >>> my_node = Node1dUpgraded(node_ini)
        >>> my_node.elements_voisins = [elem_droite, elem_gauche]
        >>> my_node.calculer_masse_wilkins()
        >>> my_node.calculer_nouvo_force()
        >>> my_node.position_relative = -1
        >>> my_node.calculer_nouvo_vitesse(1.0e-01)
        >>> print my_node.upundemi_classique
        [ 300.]
        >>> print my_node.upundemi_enrichi
        [-700.]
        >>> print my_node.upundemi
        [ 1000.]
        """
        self._upundemi_enrichi = \
            self.force_enrichi / self.masse * delta_t + self.umundemi_enrichi
        self._upundemi_classique = \
            self.force_classique / self.masse * delta_t + \
                self.umundemi_classique
        self._upundemi = \
            self.upundemi_classique + \
            self.position_relative * self.upundemi_enrichi

    def calculer_nouvo_force(self):
        """
        Calcul de la force agissant sur le noeud

        @TODO : prise en compte des CLs

        TEST UNITAIRE
        >>> import numpy as np
        >>> class element:
        ...     pass
        ...
        >>> elem_gauche = element()
        >>> elem_droite = element()
        >>> elem_gauche.coord = np.array([-0.5])
        >>> elem_droite.coord = np.array([0.5])
        >>> elem_gauche.pressure = 2.5e+09
        >>> elem_droite.pressure = 1.0e+09
        >>> node_ini = Node1d(123, section=1.0e-06)
        >>> my_node = Node1dUpgraded(node_ini)
        >>> my_node.elements_voisins = [elem_droite, elem_gauche]
        >>> for elem in my_node.elements_voisins:
        ...     print elem.coord
        [-0.5]
        [ 0.5]
        >>> my_node.calculer_nouvo_force()
        >>> print my_node.force
        None
        >>> print my_node.force_classique
        [ 1500.]
        >>> print my_node.force_enrichi
        [-3500.]
        """
        if (self.position_relative == 1):
            # Noeud à droite de la discontinuité
            pgauche = self.elements_voisins[0].pression_t_plus_dt + \
                self.elements_voisins[0].pseudo_plus_un_demi
            pgauche_enr = \
                self.elements_voisins[0]._pression_t_plus_dt_enrichi + \
                self.elements_voisins[0]._pseudo_plus_un_demi_enrichi
            pdroite = self.elements_voisins[1].pression_t_plus_dt + \
                self.elements_voisins[1].pseudo_plus_un_demi
            #
            self._force_classique[:] = (pgauche - pdroite) * self.section
            self._force_enrichi[:] = (-pgauche_enr - pdroite) * self.section
        elif (self.position_relative == -1):
            # Noeud à gauche de la discontinuité
            pgauche = self.elements_voisins[0].pression_t_plus_dt + \
                self.elements_voisins[0].pseudo_plus_un_demi
            pdroite = self.elements_voisins[1].pression_t_plus_dt + \
                self.elements_voisins[1].pseudo_plus_un_demi
            pdroite_enr = \
                self.elements_voisins[1]._pression_t_plus_dt_enrichi + \
                self.elements_voisins[1]._pseudo_plus_un_demi_enrichi
            #
            self._force_classique[:] = (pgauche - pdroite) * self.section
            self._force_enrichi[:] = (pgauche_enr - pdroite) * self.section
        self._force = None

    def incrementer(self):
        """
        Mise à jour de la vitesse et de la coordonnée du noeud
        pour passer au pas de temps suivant.
        """
        Node1d.incrementer(self)
        self._umundemi_classique[:] = self.upundemi_classique[:]
        self._umundemi_enrichi[:] = self.upundemi_enrichi[:]

if __name__ == "__main__":
    import doctest
    testres = doctest.testmod(verbose=0)
    if(testres[0] == 0):
        print "TESTS UNITAIRES : OK"
        NODE_INI = Node1d(123, section=1.0e-06)
        MY_NODE = Node1dUpgraded(NODE_INI)
        MY_NODE.infos()