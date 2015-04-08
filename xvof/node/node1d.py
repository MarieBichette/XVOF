#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un noeud en 1d
"""
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from xvof.node import Node
import numpy as np

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


class Node1d(Node):
    """
    Une classe pour les noeuds classiques dans le cas 1d
    """
    def __init__(self, indice, poz_init=np.zeros(1), vit_init=np.zeros(1),
                section=1.):
        Node.__init__(self, dim=1, index=indice, position_initiale=poz_init,
                      vitesse_initiale=vit_init)

        #
        self._section = section

    #------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    #------------------------------------------------------------

    @property
    def section(self):
        """
        Surface associée au noeud
        """
        return self._section

    @property
    def elements_voisins(self):
        """
        Liste des éléments voisins du noeud
        """
        return self._elements_voisins

    @elements_voisins.setter
    def elements_voisins(self, elems):
        """
        Setter des elements voisins. Surcharge de la méthode de Node pour
        s'assurer qu'il n y ait que deux voisins possibles et pour les trier
        de gauche à droite
        """
        if(len(elems) != 2):
            message = "En 1d seuls deux éléments peuveut être"
            message += " voisins du {}".format(self)
            raise SystemExit(message)
        self._elements_voisins = elems[:]
        self._elements_voisins =\
        sorted(self._elements_voisins, key=lambda m: m.coord[0])

    #------------------------------------------------------------
    # DEFINITIONS DES METHODES
    #------------------------------------------------------------

    def infos(self):
        """
        Affichage des informations
        """
        Node.infos(self)
        message = "==> section = {:5.4g}".format(self.section)
        print message

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
        >>> my_node = Node1d(123, section=1.0e-06)
        >>> my_node.elements_voisins = [elem_droite, elem_gauche]
        >>> for elem in my_node.elements_voisins:
        ...     print elem.coord
        [-0.5]
        [ 0.5]
        >>> my_node.calculer_nouvo_force()
        >>> print my_node.force
        [ 1500.]
        """
        self._force[:] = (self.elements_voisins[0].pressure -
            self.elements_voisins[1].pressure) * self.section

    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps supérieur

        TEST UNITAIRE
        >>> import numpy as np
        >>> class element:
        ...     pass
        ...
        >>> elem_gauche = element()
        >>> elem_gauche.coord = np.array([-0.5])
        >>> elem_gauche.masse = 3./4.
        >>> elem_gauche.pressure = 2.5e+09
        >>> elem_droite = element()
        >>> elem_droite.coord = np.array([0.5])
        >>> elem_droite.masse = 1./4.
        >>> elem_droite.pressure = 1.0e+09
        >>> my_node = Node1d(123, section=1.0e-06)
        >>> my_node.elements_voisins = [elem_droite, elem_gauche]
        >>> my_node.calculer_nouvo_force()
        >>> my_node.calculer_masse_wilkins()
        >>> my_node.calculer_nouvo_vitesse(delta_t=1.0e-01)
        >>> print my_node.upundemi
        [ 300.]
        """
        self._upundemi = self.force / self.masse * delta_t + self.umundemi


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=0)
    print "Test unitaire : OK"
    MY_NODE = Node1d(123, section=1.0e-06)
    MY_NODE.infos()