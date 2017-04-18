#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
inv operation  tests
"""

import numpy as np
from numpy.linalg import inv
import unittest

class InvMatrixTest(unittest.TestCase):
   def setUp(self) :
       self._test_matrix = np.array([[1 , 2] , [ 1 , 4 ]])
       print self._test_matrix



   def test_inverse_matrice(self):
        inverse = inv(self._test_matrix)
        inverse_solution = np.array([[2 , -1] , [ -0.5 , 0.5 ]])
        np.testing.assert_array_equal(inverse, inverse_solution)



if __name__ == '__main__':
    unittest.main()
