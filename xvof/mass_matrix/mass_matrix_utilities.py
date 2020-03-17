#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
un module pour des utilitaires numpy
"""
import numpy as np
from numpy.linalg import inv


def multiplicationMasse(matrice, vecteur):
    """
    Fonction pour faire le produit matriciel matrice * vecteur adapté pour la matrice masse sous forme de vecteur

    >>> import numpy as np
    >>> matrice  = np.array([1., 2., 3., 4.])
    >>> vecteur = np.array([1., 1./2., 1./3., 1./4.])
    >>> multiplicationMasse(matrice, vecteur)
    array([ 1.,  1.,  1.,  1.])
    >>> matrice_bis = np.array([[1., 2., 3., 4.], [2., 4., 0.5, -1], [3., 0.5, 1., -2.], [4., -1, -2., 3.]])
    >>> vecteur_bis = np.array([1., 1./2., 1./4., 1./4.])
    >>> multiplicationMasse(matrice_bis, vecteur_bis)
    array([ 3.75 ,  3.875,  3.   ,  3.75 ])
    """
    # la matrice est en fait un vecteur ligne ou colonne--> produit terme à terme
    if len(matrice.shape) == 1 or matrice.ndim == 2 and 1 in matrice.shape:
        result = matrice * vecteur
    else:  # matrice est unevraie matrice (dim >=2) --> produit matriciel
        result = np.dot(matrice, vecteur)
    return result


def inverseMasse(matrice):
    """
    MassMatrix de type MassMatrix
    Fonction pour faire inverse la matrice masse, qu'elle soit sous forme de vecteur ou de matrice
    """
    if len(matrice.shape) == 1 or (matrice.ndim == 2 and 1 in matrice.shape):
        result = np.zeros([len(matrice), 1])
        for ind_node in xrange(len(matrice)):
            result[ind_node] = 1. / matrice[ind_node]
    else:  # la matrice est une vraie matrice (dim >=2)
        result = inv(matrice)
    return result


def lump_matrix(matrix):
    """
    Condense la matrice de masse avec la méthode de Menouillard
    (on somme sur toute la ligne pour obtenir une matrice diagonale)
    """
    lumped_matrix = np.zeros((matrix.shape[0], matrix.shape[1]))
    for iligne in range(matrix.shape[0]):
        lumped_matrix[iligne, iligne] = sum(matrix[iligne, :])
    return lumped_matrix


class SymNDArray(np.ndarray):
    """
    Une classe pour les matrices symétriques
    """
    def __setitem__(self, (i, j), value):
        super(SymNDArray, self).__setitem__((i, j), value)
        super(SymNDArray, self).__setitem__((j, i), value)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
