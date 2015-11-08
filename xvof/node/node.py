#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module définissant la classe Node
"""

from abc import abstractmethod
import numpy as np

# @Todo : Créer un set (variable de classe) qui contient l'ensemble
# des indices de tous les noeuds pour éviter les doublons


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
    def __init__(self, dim=1, nbr_of_nodes, position_initiale=None, vitesse_initiale=None):
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
	# Le vecteur est un vecteur dont les lignes sont les noeuds et les colonnes les coordonées selon
	# les différentes dimensions
	self.__shape = [nbr_of_nodes, dim]
        # Les autres attributs ne sont pas publics mais restent accessibles et
        # modifiables par les classes filles
        if position_initiale is None:
            position_initiale = np.zeros(self.__shape, dtype=np.float64, order='C')
        elif np.shape(position_initiale) != self.__shape:
            message = "Node() : La dimension ({}) du vecteur position_initiale "\
                .format(np.shape(position_initiale))
            message += "est incorrecte!"
            raise SystemExit(message)
        if vitesse_initiale is None:
            vitesse_initiale = np.zeros(self.__shape, dtype=np.float64, order='C')
        elif np.shape(vitesse_initiale) != self.__shape:
            message = "Node() : La dimension ({}) du vecteur vitesse_initiale "\
                .format(np.shape(vitesse_initiale))
            message += "est incorrecte!"
            raise SystemExit(message)
        #
        self._xt = np.array(position_initiale)
        self._umundemi = np.array(vitesse_initiale)
        self._xtpdt = np.zeros(self.__shape, dtype=np.float64, order='C')
        self._upundemi = np.array(vitesse_initiale)
        self._masse = np.zeros([self.__shape[0], 1], dtype=np.float64, order='C')
        self._invmasse = np.zeros([self.__shape[0], 1], dtype=np.float64, order='C')
        self._force = np.zeros([self.__shape[0], 1], dtype=np.float64, order='C')

    property
    def coordt(self):
        """
        Positions des noeuds au temps t

	:return: positions des noeuds au temps t
	:rtype: numpy.array([nbr_of_nodes, dim], dtype=np.float64, order='C')
        """
        return self._xt

    @property
    def coordtpdt(self):
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
    def invmasse(self):
        """
        Inverse des masses nodales

	:return: vecteur des inverses des masses nodales
	:rtype: numpy.array([nbr_of_nodes, 1], dtype=np.float64, order='C')
        """
        return self._invmasse

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

    @property:
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
        message += "==> coordonnées à t = {}\n".format(self.coordt[index])
        message += "==> coordonnées à t+dt = {}\n".format(self.coordtpdt[index])
        message += "==> vitesse à t-1/2 = {}\n".format(self.umundemi[index])
        message += "==> vitesse à t+1/2 = {}\n".format(self.upundemi[index])
        message += "==> masse = {:5.4g}\n".format(self.masse[index])
        message += "==> force = {}".format(self.force[index])
        print message

    def calculer_masse_wilkins(self, elements_voisins):
        """
        Calcule la masse associée au noeud par moyenne arithmétique de la
        masse des éléments voisins (méthode Wilkins)
        """
        for elem in elements_voisins:
            self._masse += elem.masse / elem.nbr_noeuds
            self._invmasse = 1. / self._masse

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
    def calculer_nouvo_force(self, *args, **kwargs):
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
