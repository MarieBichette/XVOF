#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining the mass matrix 1D
"""

import numpy as np

from xvof.mass_matrix.mass_matrix import MassMatrix


class OneDimensionMassMatrix(MassMatrix):
    """
    A class for mass matrix 1d
    """

    def __init__(self, number_of_nodes):
        MassMatrix.__init__(self, number_of_nodes)
        self._inv_mass_matrix = np.zeros([self.number_of_nodes, 1], dtype=np.float64, order='C')

    def compute_mass_matrix(self, topologie, vecteur_masse_elements, vecteur_nb_noeuds_par_element):
        """
        Compute the mass matrix (methode de Wilkins)
        self.mass_matix is type array 1D (= vecteur)
        """
        self._mass_matrix = MassMatrix.compute_mass_matrix(self, topologie, vecteur_masse_elements,
                                                           vecteur_nb_noeuds_par_element)

    @property
    def mass_matrix_value(self):
        """Accessor on _mass_matrix"""
        return self._mass_matrix


