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

        self.value = 10.
        self.criterion = NonLocalStressCriterion(self.value)

        self.cells = OneDimensionHansboEnrichedCell(10)


    def test_check_criterion_rupture_simple(self):
        """
        Test of the method check_criterion of MinimumPressureCriterion
        """
        self.cells._stress = np.array([[1, ], [4, ], [25, ], [7, ], [-2, ],
                                       [10, ], [26, ], [1, ], [0, ], [15, ]])
        result = self.criterion.check_criterion(self.cells)
        exact = np.array([False, False, True, False, False,
                          True, True, False, False, False])
        np.testing.assert_allclose(exact, result)


if __name__ == '__main__':
    unittest.main()
