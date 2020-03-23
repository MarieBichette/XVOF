# -*- coding: utf-8 -*-
"""
inv operation tests

:todo: this unittest seems useless
"""
import unittest
import numpy as np
from numpy.linalg import inv

class InvMatrixTest(unittest.TestCase):
    """This class tests the matrix inversion process"""
    def setUp(self):
        self._test_matrix = np.array([[1, 2], [1, 4]])
        print(self._test_matrix)

    def test_inverse_matrice(self):
        """Test the inverse of a matrix"""
        inverse = inv(self._test_matrix)
        inverse_solution = np.array([[2, -1], [-0.5, 0.5]])
        np.testing.assert_array_equal(inverse, inverse_solution)


if __name__ == '__main__':
    unittest.main()
