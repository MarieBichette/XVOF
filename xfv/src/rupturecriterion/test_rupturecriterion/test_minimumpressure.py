# -*- coding: utf-8 -*-
"""
Classe de test du module node
"""
import unittest
import unittest.mock as mock
import numpy as np

from xfv.src.cell.cell import Cell
from xfv.src.rupturecriterion.minimumpressure import MinimumPressureCriterion


class MinimumPressureTest(unittest.TestCase):
    """
    Test case utilis� pour test les fonctions du module 'ConstantPressure'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        self._prupturemin = 50.
        self._rupturecriterion = MinimumPressureCriterion(self._prupturemin)

        self.cell = mock.MagicMock(Cell)
        self.cell.pressure = mock.PropertyMock()
        self.cell.pressure.new_value = np.zeros([1000])
        self.cell.pressure.new_value[500] = 25.
        self.cell.pressure.new_value[501] = 75.

    def test_check_criterion_rupture(self):
        """
        Test of the method check_criterion of MinimumPressureCriterion
        """
        retour = self._rupturecriterion.check_criterion(self.cell)
        self.assertFalse(retour[501])
        self.assertTrue(retour[500])


if __name__ == '__main__':
    unittest.main()
