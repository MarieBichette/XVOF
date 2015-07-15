#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant un noeud
"""

from abc import abstractmethod
import numpy as np

# @Todo : Créer un set (variable de classe) qui contient l'ensemble
# des indices de tous les noeuds pour éviter les doublons


class Node(object):
    """
    Une classe pour les noeuds
    """
    # pylint: disable-msg=R0902
    # 10 attributs : cela semble raisonnable pour ce cas
    def __init__(self, dim=1, index=-1,
                 position_initiale=None,
                 vitesse_initiale=None):
        # __dimension doit rester privé même pour les classes filles
        # un noeud consacré aux simulations 1d ne peut pas changer sa dimension
        self.__dimension = dim
        # Les autres attributs ne sont pas publics mais restent accessibles et
        # modifiables par les classes filles
        if (not isinstance(index, int)):
            raise TypeError("L'indice du noeud doit être un entier!")
        self._index = index
        #
        if (position_initiale is None):
            position_initiale = np.zeros(self.__dimension, dtype=float)
        elif (np.shape(position_initiale) != (self.__dimension,)):
            message = "Node() : La dimension ({}) du vecteur position_initiale "\
                .format(np.shape(position_initiale))
            message += "est incorrecte!"
            raise SystemExit(message)
        if (vitesse_initiale is None):
            vitesse_initiale = np.zeros(self.__dimension, dtype=float)
        elif (np.shape(vitesse_initiale) != (self.__dimension,)):
            message = "La dimension du vecteur position_initiale "
            message += "est incorrecte!"
            raise SystemExit(message)
        #
        self._xt = np.array(position_initiale)
        self._umundemi = np.array(vitesse_initiale)
        self._xtpdt = np.zeros(self.__dimension, dtype=float)
        self._upundemi = np.array(vitesse_initiale)
        self._masse = 0.
        self._force = np.zeros(self.__dimension, dtype=float)

    @property
    def index(self):
        """
        Indice global du noeud
        """
        return self._index

    @property
    def coordt(self):
        """
        Position du noeud au temps t
        """
        return self._xt

    @property
    def coordtpdt(self):
        """
        Position du noeud au temps t + dt
        """
        return self._xtpdt

    @property
    def umundemi(self):
        """
        Vitesse au demi pas de temps précédent
        """
        return self._umundemi

    @property
    def upundemi(self):
        """
        Vitesse au demi pas de temps suivant
        """
        return self._upundemi

    @property
    def masse(self):
        """
        Masse nodale
        """
        return self._masse

    @property
    def force(self):
        """
        Force nodale
        """
        return self._force

    @property
    def dimension(self):
        """
        Dimension associée
        """
        return self.__dimension

    def __str__(self):
        message = "NOEUD {:4d} ".format(self.index)
        message += "(dimension : {:1d})".format(self.__dimension)
        return message

    def infos(self):
        """
        Affichage des informations concernant le noeud
        """
        message = "{} {:4d}\n".format(self.__class__, self.index)
        message += "==> coordonnées à t = {}\n".format(self.coordt)
        message += "==> coordonnées à t+dt = {}\n".format(self.coordtpdt)
        message += "==> vitesse à t-1/2 = {}\n".format(self.umundemi)
        message += "==> vitesse à t+1/2 = {}\n".format(self.upundemi)
        message += "==> masse = {:5.4g}\n".format(self.masse)
        message += "==> force = {}".format(self.force)
        print message

    def calculer_masse_wilkins(self, elements_voisins):
        """
        Calcule la masse associée au noeud par moyenne arithmétique de la
        masse des éléments voisins (méthode Wilkins)
        """
        for elem in elements_voisins:
            self._masse += elem.masse / elem.nbr_noeuds

    def calculer_nouvo_coord(self, delta_t=1.0):
        """
        Calcul de la coordonnée au temps t+dt

        @param delta_t : pas de temps
        """
        self._xtpdt = self.coordt + self.upundemi * delta_t

    def incrementer(self):
        """
        Mise à jour de la vitesse et de la coordonnée du noeud
        pour passer au pas de temps suivant.
        """
        self._umundemi[:] = self.upundemi[:]
        self._xt[:] = self.coordtpdt[:]

    @abstractmethod
    def calculer_nouvo_force(self):
        """
        Calcul de la force agissant sur le noeud
        """

    @abstractmethod
    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps supérieur
        """

if __name__ == "__main__":
    print "Ceci est uniquement un module!"
