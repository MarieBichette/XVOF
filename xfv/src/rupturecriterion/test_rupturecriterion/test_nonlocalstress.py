# -*- coding: utf-8 -*-
"""
Classe de test du module NonLocalCriterion
"""
import unittest
import unittest.mock as mock
import os
import numpy as np

from xfv.src.cell.one_dimension_enriched_cell_hansbo import OneDimensionHansboEnrichedCell
from xfv.src.cell.cell import Cell
from xfv.src.rupturecriterion.nonlocalstress import NonLocalStressCriterion
from xfv.src.rupturecriterion.nonlocalstressweight import ArithmeticWeight
from xfv.src.data.data_container import DataContainer


class NonLocalCriterionTest(unittest.TestCase):
    """
    Test case utilis� pour test les fonctions du module 'NonLocalCriterion'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_enrichment_epp.json")
        DataContainer(data_file_path)

    def test_check_criterion_rupture_1(self):
        """
        Test of the method check_criterion of NonLocalStressCriterion
        WITHOUT ENRICHMENT
        """
        # configuration d'un mock cell
        config = {'cell_in_target': np.ones(10, dtype=bool),
                  'enriched': np.zeros(10, dtype=bool),
                  'coordinates_x': np.array([[0., ], [1., ], [2., ], [3., ], [4., ],
                                             [5., ], [6., ], [7., ], [8., ], [9., ]]),
                  'stress_xx': np.array([[1, ], [4, ], [25, ], [7, ], [-2, ], [10, ], [26, ], [1, ], [0, ], [15, ]]),
                  'enr_coordinates_x': np.zeros([10, 1]),
                  'enr_stress_xx': np.zeros([10, 1]),
                  'toto': 12}
        patcher = mock.patch('xfv.src.cell.one_dimension_enriched_cell_hansbo.OneDimensionHansboEnrichedCell',
                             spec=OneDimensionHansboEnrichedCell, **config)
        mock_cells = patcher.start()

        self.criterion = NonLocalStressCriterion(10., ArithmeticWeight(1.1))
        result = self.criterion.check_criterion(mock_cells)
        exact = np.array([False, True, True, True, False, True, True, False, False, False])
        np.testing.assert_allclose(exact, result)

    def test_check_criterion_rupture_2(self):
        """
        Test of the method check_criterion of NonLocalStressCriterion
        WITHOUT ENRICHMENT
        """
        # configuration d'un mock cell
        config = {'cell_in_target': np.ones(10, dtype=bool),
                  'enriched': np.zeros(10, dtype=bool),
                  'coordinates_x': np.array([[0., ], [1., ], [2., ], [3., ], [4., ],
                                             [5., ], [6., ], [7., ], [8., ], [9., ]]),
                  'stress_xx': np.array([[1, ], [4, ], [25, ], [7, ], [-2, ], [10, ], [26, ], [1, ], [0, ], [15, ]]),
                  'enr_coordinates_x': np.zeros([10, 1]),
                  'enr_stress_xx': np.zeros([10, 1]),
                  'toto': 12}
        patcher = mock.patch('xfv.src.cell.one_dimension_enriched_cell_hansbo.OneDimensionHansboEnrichedCell',
                             spec=OneDimensionHansboEnrichedCell, **config)
        mock_cells = patcher.start()

        self.criterion = NonLocalStressCriterion(10., ArithmeticWeight(2.1))
        result = self.criterion.check_criterion(mock_cells)
        exact = np.array([True, False, False, False, True, False, False, True, True, False])
        np.testing.assert_allclose(exact, result)


if __name__ == '__main__':
    unittest.main()
