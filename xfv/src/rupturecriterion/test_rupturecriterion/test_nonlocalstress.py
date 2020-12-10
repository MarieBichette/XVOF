# -*- coding: utf-8 -*-
"""
Classe de test du module node
"""
import unittest
import unittest.mock as mock
import os
import numpy as np

from xfv.src.cell.one_dimension_enriched_cell_hansbo import OneDimensionHansboEnrichedCell
from xfv.src.rupturecriterion.nonlocalstress import NonLocalStressCriterion
from xfv.src.data.data_container import DataContainer


class MinimumPressureTest(unittest.TestCase):
    """
    Test case utilis� pour test les fonctions du module 'ConstantPressure'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_enrichment_epp.json")
        DataContainer(data_file_path)

        self.cells = OneDimensionHansboEnrichedCell(10)


    def test_check_criterion_rupture_1(self):
        """
        Test of the method check_criterion of MinimumPressureCriterion
        """
        self.cells._coordinates_x = np.array([[0., ], [1., ], [2., ], [3., ], [4., ],
                                              [5., ], [6., ], [7., ], [8.,], [9., ]])
        self.cells._stress = np.array([[1, 0, 0], [4, 0, 0], [25, 0, 0], [7, 0, 0], [-2, 0, 0],
                                       [10, 0, 0], [26, 0, 0], [1, 0, 0], [0, 0, 0], [15, 0, 0]])
        self.criterion = NonLocalStressCriterion(10., 1.)
        result = self.criterion.check_criterion(self.cells)
        exact = np.array([False, False, True, False, False, True, True, False, False, True])
        np.testing.assert_allclose(exact, result)

    def test_check_criterion_rupture_2(self):
        """
        Test of the method check_criterion of MinimumPressureCriterion
        """
        self.cells._coordinates_x = np.array([[0., ], [1., ], [2., ], [3., ], [4., ],
                                              [5., ], [6., ], [7., ], [8.,], [9., ]])
        self.cells._stress = np.array([[1, 0, 0], [4, 0, 0], [25, 0, 0], [7, 0, 0], [-2, 0, 0],
                                       [10, 0, 0], [26, 0, 0], [1, 0, 0], [0, 0, 0], [15, 0, 0]])
        self.criterion = NonLocalStressCriterion(10., 2.)
        result = self.criterion.check_criterion(self.cells)
        exact = np.array([False, True, True, True, False, True, True, False, False, False])
        np.testing.assert_allclose(exact, result)


if __name__ == '__main__':
    unittest.main()
