#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module ImposedPressure
"""
import unittest
import xvof.src.rupturetreatment.imposedpressure as IP
import numpy as np

from xvof.src.cell.one_dimension_cell import OneDimensionCell


class ImposedPressureTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'Imposed Pressure'
    Marche bien mais appelle des méthodes de Cell et OneDimensionCell non vérifiées
    """

    def setUp(self):
        """
        Préparation du test
        """
        self._pressure = 0.

        self.my_imposed_pressure = IP.ImposedPressure(self._pressure)

        self.cells = OneDimensionCell(1000)
        self.cells.pressure_field[:] = 1.

        self.ruptured_cells = np.ndarray(1000, dtype=np.bool, order = 'C')
        self.ruptured_cells[:] = False
        self.ruptured_cells[500] = True

    def test_applyTreatment(self):
        """
        Teste la méthode applyTreatment for ImposedPressure
        """
        self.my_imposed_pressure.applyTreatment(self.cells, self.ruptured_cells)
        self.cells.increment_variables()
        self.assertAlmostEqual(self.cells.pressure_field[500], self._pressure)


if __name__ == '__main__':
    unittest.main()
