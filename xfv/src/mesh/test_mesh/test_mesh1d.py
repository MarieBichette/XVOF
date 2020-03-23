# -*- coding: utf-8 -*-
"""
Classe de test du module Mesh1d
"""
import numpy as np
import unittest
import unittest.mock as mock

from xfv.src.mesh.mesh1d import Mesh1d

class Mesh1dTest(unittest.TestCase):
    """
    Test case utilis� pour test les fonctions du module 'Mesh1d'
    """
    def setUp(self):
        """
        Pr�paration des tests unitaires
        """
        coord_init = np.array([0., 1., 2., 3.])
        vit_init = np.zeros(4)
        self.my_mesh = Mesh1d(coord_init, vit_init)








if __name__ == '__main__':
    unittest.main()