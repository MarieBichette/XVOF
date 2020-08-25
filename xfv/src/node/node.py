#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module defining Node class
"""

from abc import abstractmethod
import numpy as np


class Node:
    """
    Un objet Node repr�sente l'ensemble des noeuds du maillages. Ses diff�rents membres sont
    essentiellement des vecteurs de nbr_of_nodes lignes. Plusieurs colonnes peuvent �tre
    pr�sentes selon la dimension du probl�me � traiter.

    L'organisation en m�moire est comme en C/C++ c'est � dire 'row wise'. C'est pour cette raison
    que les lignes de chacun des vecteurs repr�sentent les noeuds. Ce faisant on a par exemple
    tous les X des noeuds contigus en m�moire. (vectorisation, localisation spatiale)
    """

    # pylint: disable-msg=R0902
    # 10 attributs : cela semble raisonnable pour ce cas
    def __init__(self, nbr_of_nodes, position_initiale, dim=1, vitesse_initiale=None):
        """
        :param dim: dimension du probl�me � traiter (par d�faut 1)
        :param nbr_of_nodes: nombre de noeuds du probl�me
        :param position_initiale: vecteur des positions initiales
        :param vitesse_initiale: vecteur des vitesses initiales
        :type dim: int
        :type nbr_of_nodes: int
        :type position_initiale: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        :type vitesse_initiale: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        """
        # Le vecteur est un vecteur dont les lignes sont les noeuds et les colonnes
        # les coordon�es selon les diff�rentes dimensions
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
    def xt(self):  # pylint: disable=invalid-name
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
        Vitesses au demi pas de temps pr�c�dent

        :return: vitesses des noeuds au demi pas de temps pr�c�dent
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
        Dimension du probl�me

        :return: dimension du probl�me
        :rtype: int
        """
        return self.__shape[1]

    @property
    def number_of_nodes(self):
        """
        Nombre de noeuds du probl�me

        :return: dimension du probl�me
        :rtype: int
        """
        return self.__shape[0]

    def __str__(self):
        message = "Nombre de noeuds {:8d} ".format(self.__shape[0])
        message += "Dimension du probl�me : {:1d})".format(self.__shape[1])
        return message

    def infos(self, index):
        """
        Affichage des informations concernant le noeud d'indice index

        :param index: indice du noeud � afficher
        :type index: int
        """
        message = "{} {:4d}\n".format(self.__class__, index)
        message += "==> coordonn�es � t = {}\n".format(self.xt[index])
        message += "==> coordonn�es � t+dt = {}\n".format(self.xtpdt[index])
        message += "==> vitesse � t-1/2 = {}\n".format(self.umundemi[index])
        message += "==> vitesse � t+1/2 = {}\n".format(self.upundemi[index])
        message += "==> masse = {:5.4g}\n".format(self.masse[index])
        message += "==> force = {}".format(self.force[index])
        print(message)

    def compute_new_coodinates(self, mask: np.array, delta_t: float):
        """
        Calcul de la coordonn�e au temps t+dt
        :param mask : mask to select some specific nodes
        :param delta_t: time step
        """
        self._xtpdt[mask] = self.xt[mask] + self.upundemi[mask] * delta_t

    def increment(self):
        """
        Mise � jour de la vitesse et de la coordonn�e des noeuds
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
        Calcul de la vitesse au demi pas de temps sup�rieur
        """


if __name__ == "__main__":
    print("Ceci est uniquement un module!")
