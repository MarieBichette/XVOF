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
    def __init__(self, indice, poz_init=np.zeros(1), vit_init=np.zeros(1),
                section=1.):
        Node1d.__init__(self, indice, poz_init=poz_init,
                      vit_init=vit_init, section=section)

        self._upundemi = vit_init[:]
        self._force = np.zeros(1, dtype=float)
        #
        self._umundemi_classique = vit_init[:]
        self._upundemi_classique = vit_init[:]
        self._force_classique = np.zeros(1, dtype=float)
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

    def initialize(self, vitesse_t_m, vitesse_t_p, force):
        """
        Initialisation des champs classiques

        @param vitesse_t_m : vecteur vitesse au demi pas de temps précédent
        @param vitesse_t_p : vecteur vitesse au demi pas de temps suivant
        @param force : vecteur force

        TEST UNITAIRE
        >>> import numpy as np
        >>> MY_NODE = Node1dUpgraded(123, section=1.0e-06)
        >>> MY_NODE.initialize([-1.0], [2.5], [3.0e+04])
        >>> print MY_NODE.umundemi_classique
        [-1.]
        >>> print MY_NODE.upundemi_classique
        [ 2.5]
        >>> print MY_NODE.force_classique
        [ 30000.]
        >>> MY_NODE2 = Node1dUpgraded(123, section=1.0e-06)
        >>> MY_NODE2.initialize(np.array([-1.0]), np.array([2.5]),\
                                np.array([3.0e+04]))
        >>> print MY_NODE2.umundemi_classique
        [-1.]
        >>> print MY_NODE2.upundemi_classique
        [ 2.5]
        >>> print MY_NODE2.force_classique
        [ 30000.]
        """
        self._umundemi_classique = np.array(vitesse_t_m[:])
        self._upundemi_classique = np.array(vitesse_t_p[:])
        self._force_classique = np.array(force[:])

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
        >>> my_node = Node1dUpgraded(123, section=1.0e-06)
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
        >>> my_node = Node1dUpgraded(123, section=1.0e-06)
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
        self._force_classique[:] = (self.elements_voisins[0].pressure -
            self.elements_voisins[1].pressure) * self.section
        self._force_enrichi[:] = (-self.elements_voisins[0].pressure -
            self.elements_voisins[1].pressure) * self.section
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
        MY_NODE = Node1dUpgraded(123, section=1.0e-06)
        MY_NODE.infos()