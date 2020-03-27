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
        message = ""
        # Test dimension of array 2d
        try:
            CohesiveLaw(np.array([1, 2, 3]))
        except AssertionError as e:
            message = e.args[0]
        self.assertEqual(message, "array should be 2D")
        # Test size of the array
        try:
            CohesiveLaw(np.array([[1, 2, 3], [4, 5, 6]]))
        except AssertionError as e:
            message = e.args[0]
        self.assertEqual(message, "array should be size (x, 2)")

        # Test value of first value of separation = 0
        try:
            CohesiveLaw(np.array([[4, 2], [5, 0]]))
        except AssertionError as e:
            message = e.args[0]
        self.assertEqual(message, "first value of separation should be 0.")

        # Test value of stress at critical separation = 0
        try:
            CohesiveLaw(np.array([[0, 2], [1, 1]]))
        except AssertionError as e:
            message = e.args[0]
        self.assertEqual(message, "last value of stress should be 0.")

        # Test separation are croissant in array
        try:
            CohesiveLaw(np.array([[0, 2], [2, 1], [1, 0]]))
        except AssertionError as e:
            message = e.args[0]
        self.assertEqual(message, "separation is not sorted")

    def test_compute_cohesive_force(self):
        """
        Test of the method compute_cohesive_force of module CohesiveLaw
        """
        linear_law = CohesiveLaw(np.array([[0., 10.], [5., 0.]]))
        # Test linear law
        result = linear_law.compute_cohesive_force(4.)
        expected = 2.
        self.assertEqual(result, expected)

        bilinear_law = CohesiveLaw(np.array([[0., 10.], [3., 10.], [5., 0.]]))
        # Test bilinear law, part 1
        result = bilinear_law.compute_cohesive_force(2.)
        expected = 10.
        self.assertEqual(result, expected)

        # Test bilinear law, part 2
        result = bilinear_law.compute_cohesive_force(4.)
        expected = 5.
        self.assertEqual(result, expected)

        trilinear_law = CohesiveLaw(np.array([[0., 10.], [1., 8.], [3., 8], [5., 0.]]))
        # Test trilinear law, part 1
        result = trilinear_law.compute_cohesive_force(0.5)
        expected = 9.
        self.assertEqual(result, expected)

        # Test trilinear law, part 2
        result = trilinear_law.compute_cohesive_force(2.)
        expected = 8.
        self.assertEqual(result, expected)

         # Test trilinear law, part 3
        result = trilinear_law.compute_cohesive_force(4.)
        expected = 4.
        self.assertEqual(result, expected)




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
if __name__ == '__main__':
    unittest.main()
