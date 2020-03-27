#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Classe de test du module cohesive law
"""
import unittest
import os
import numpy as np

from xfv.src.cohesive_model.cohesive_law import CohesiveLaw
from xfv.src.data.data_container import DataContainer


class CohesiveLawTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'CohesiveLaw' for linear cohesive law
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        pass

    def tearDown(self):
        """
        Operations to be done after completing all the tests in the class
        """
        pass

    def test_init(self):
        """
        Test of the creation of CohesiveLaw
        :return:
        """
        # Test dimension of array

        # Test value of first value of separation = 0

        # Test value of stress at critical separation = 0

        # Test separation are croissant in array

    def test_compute_cohesive_force(self):
        """
        Test of the method compute_cohesive_force of module CohesiveLaw
        :return:
        """
        CohesiveLaw(np.array([[0., 10.], [5., 0.]]))




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
if __name__ == '__main__':
    unittest.main()
