#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module defining Node class
"""

from abc import abstractmethod
import numpy as np


class Node(object):
    """
    Un objet Node représente l'ensemble des noeuds du maillages. Ses différents membres sont
    essentiellement des vecteurs de nbr_of_nodes lignes. Plusieurs colonnes peuvent être
    présentes selon la dimension du problème à traiter.

    L'organisation en mémoire est comme en C/C++ c'est à dire 'row wise'. C'est pour cette raison
    que les lignes de chacun des vecteurs représentent les noeuds. Ce faisant on a par exemple
    tous les X des noeuds contigus en mêmoire. (vectorisation, localisation spatiale)
    """

    # pylint: disable-msg=R0902
    # 10 attributs : cela semble raisonnable pour ce cas
    def __init__(self, nbr_of_nodes, position_initiale, dim=1, vitesse_initiale=None):
        """
        :param dim: dimension du problème à traiter (par défaut 1)
        :param nbr_of_nodes: nombre de noeuds du problème
        :param position_initiale: vecteur des positions initiales
        :param vitesse_initiale: vecteur des vitesses initiales
        :type dim: int
        :type nbr_of_nodes: int
        :type position_initiale: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        :type vitesse_initiale: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        """
        # Le vecteur est un vecteur dont les lignes sont les noeuds et les colonnes
        # les coordonées selon les différentes dimensions
        self.__shape = (nbr_of_nodes, dim)
        # Les autres attributs ne sont pas publics mais restent accessibles et
        # modifiables par les classes filles
        if position_initiale.shape != self.__shape:
            message = "Node() : La dimension ({}) du vecteur position_initiale ".format(
                np.shape(position_initiale))
            message += "est incorrecte (!= {})!".format(self.__shape)
            raise SystemExit(message)
        if vitesse_initiale is None:
            vitesse_initiale = np.zeros(self.__shape, dtype=np.float64, order='C')
        elif np.shape(vitesse_initiale) != self.__shape:
            message = "Node() : La dimension ({}) du vecteur vitesse_initiale " \
                .format(np.shape(vitesse_initiale))
            message += "est incorrecte (!= {})!".format(self.__shape)
            raise SystemExit(message)
        #
        self._xt = np.array(position_initiale)
        self._umundemi = np.array(vitesse_initiale)
        self._xtpdt = np.zeros(self.__shape, dtype=np.float64, order='C')
        self._upundemi = np.array(vitesse_initiale)

        self._masse = np.zeros(self.__shape, dtype=np.float64, order='C')

        self._force = np.zeros([self.__shape[0], dim], dtype=np.float64, order='C')

        self._enriched = np.empty(nbr_of_nodes)
        self._enriched[:] = False

    @property
    def enriched(self):
        """
        Returns an array of the status of all nodes
        False = classical node
        True = enriched node
        :return:  numpy.array([nbr_of_nodes], dtype=bool)
        """
        return self._enriched

    @property
    def xt(self):
        """
        Positions des noeuds au temps t

        :return: positions des noeuds au temps t
        :rtype: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        """
        return self._xt

    @property
    def xtpdt(self):
        """
        Positions des noeuds au temps t + dt

        :return: positions des noeuds au temps t + dt
        :rtype: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        """
        return self._xtpdt

    @property
    def umundemi(self):
        """
        Vitesses au demi pas de temps précédent

        :return: vitesses des noeuds au demi pas de temps précédent
        :rtype: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        """
        return self._umundemi

    @property
    def upundemi(self):
        """
        Vitesses au demi pas de temps suivant

        :return: vitesses des noeuds au demi pas de temps suivant
        :rtype: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        """
        return self._upundemi

    @property
    def masse(self):
        """
        Masses nodales

        :return: vecteur des masses nodales
        :rtype: numpy.array([nbr_of_nodes, 1], dtype=np.float64, order='C')
        """
        return self._masse

    @property
    def force(self):
        """
        Forces nodales

        :return: vecteur des forces nodales
        :rtype: numpy.array([nbr_of_nodes, 1], dtype=np.float64, order='C')
        """
        return self._force

    @property
    def dimension(self):
        """
        Dimension du problème

        :return: dimension du problème
        :rtype: int
        """
        return self.__shape[1]

    @property
    def number_of_nodes(self):
        """
        Nombre de noeuds du problème

        :return: dimension du problème
        :rtype: int
        """
        return self.__shape[0]

    def __str__(self):
        message = "Nombre de noeuds {:8d} ".format(self.__shape[0])
        message += "Dimension du problème : {:1d})".format(self.__shape[1])
        return message

    def infos(self, index):
        """
        Affichage des informations concernant le noeud d'indice index

        :param index: indice du noeud à afficher
        :type index: int
        """
        message = "{} {:4d}\n".format(self.__class__, index)
        message += "==> coordonnées à t = {}\n".format(self.xt[index])
        message += "==> coordonnées à t+dt = {}\n".format(self.xtpdt[index])
        message += "==> vitesse à t-1/2 = {}\n".format(self.umundemi[index])
        message += "==> vitesse à t+1/2 = {}\n".format(self.upundemi[index])
        message += "==> masse = {:5.4g}\n".format(self.masse[index])
        message += "==> force = {}".format(self.force[index])
        print message

    def compute_new_coodinates(self, delta_t):
        """
        Calcul de la coordonnée au temps t+dt

        :param delta_t: pas de temps
        :type delta_t: float
        """
        self._xtpdt = self.xt + self.upundemi * delta_t

    def increment(self):
        """
        Mise à jour de la vitesse et de la coordonnée des noeuds
        pour passer au pas de temps suivant.
        """
        self._umundemi[:] = self.upundemi[:]
        self._xt[:] = self.xtpdt[:]

    @abstractmethod
    def compute_new_force(self, *args, **kwargs):
        """
        Calcul de la force agissant sur le noeud
        """

    @abstractmethod
    def compute_new_velocity(self, delta_t, mask, matrice_masse):
        """
        Calcul de la vitesse au demi pas de temps supérieur
        """


if __name__ == "__main__":
    print "Ceci est uniquement un module!"
