#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module node
"""
import numpy as np
import unittest
import mock
from xvof.cell.cell import Cell
from xvof.rupturecriterion.minimumpressure import MinimumPressureCriterion

class MinimumPressureTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'ConstantPressure'
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

    def test_checkCriterion_rupture(self):
        retour = self._rupturecriterion.checkCriterion(self.cell)
        self.assertFalse(retour[501])
        self.assertTrue(retour[500])


if __name__ == '__main__':
    unittest.main()