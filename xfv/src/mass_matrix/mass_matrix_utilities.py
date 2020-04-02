#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
un module pour des utilitaires numpy
"""
import numpy as np
from numpy.linalg import inv


def multiplication_masse(matrix, vector):
    """
    Fonction pour faire le produit matriciel matrice * vecteur adapté pour la matrice masse sous 
    forme de vecteur
    :param matrix : matrix (array multiD)
    :param vector: vector (array 1D)

    >>> import numpy as np
    >>> matrice  = np.array([1., 2., 3., 4.])
    >>> vecteur = np.array([1., 1./2., 1./3., 1./4.])
    >>> multiplication_masse(matrice, vecteur)
    array([ 1.,  1.,  1.,  1.])
    >>> matrice_bis = \
    np.array([[1., 2., 3., 4.], [2., 4., 0.5, -1], [3., 0.5, 1., -2.], [4., -1, -2., 3.]])
    >>> vecteur_bis = np.array([1., 1./2., 1./4., 1./4.])
    >>> multiplication_masse(matrice_bis, vecteur_bis)
    array([ 3.75 ,  3.875,  3.   ,  3.75 ])
    """
    # matrix is a vector (line or column) => term by term product
    if len(matrix.shape) == 1 or matrix.ndim == 2 and 1 in matrix.shape:
        result = matrix * vector
    else:  # matrix is a true matrix (dim >=2) --> matrix vector product
        result = np.dot(matrix, vector)
    return result


def inverse_masse(matrix):
    """
    MassMatrix de type MassMatrix
    Fonction pour faire inverse la matrice masse, qu'elle soit sous forme de vecteur ou de matrice
    :param matrix: matrix to inverse
    """
    if len(matrix.shape) == 1 or (matrix.ndim == 2 and 1 in matrix.shape):
        result = np.zeros([len(matrix), 1])
        for ind_node in range(len(matrix)):
            result[ind_node] = 1. / matrix[ind_node]
    else:  # matrix is a true matrix (dim >=2)
        result = inv(matrix)
    return result


def lump_matrix(matrix):
    """
    Condense la matrice de masse avec la méthode de Menouillard
    (on somme sur toute la ligne pour obtenir une matrice diagonale)
    :param matrix: matrix to lump
    """
    lumped_matrix = np.zeros((matrix.shape[0], matrix.shape[1]))
    for iligne in range(matrix.shape[0]):
        lumped_matrix[iligne, iligne] = sum(matrix[iligne, :])
    return lumped_matrix


class SymNDArray(np.ndarray):
    """
    Une classe pour les matrices symétriques
    """
    def __setitem__(self, xxx_todo_changeme, value):
        """
        Build the symmetric matrix
        :param xxx_todo_changeme:
        :param value:
        :return:
        """
        (i, j) = xxx_todo_changeme
        super(SymNDArray, self).__setitem__((i, j), value)
        super(SymNDArray, self).__setitem__((j, i), value)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
