#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Classe de test du module ImposedPressure
"""
import unittest
import os
import numpy as np

from xfv.src.rupturetreatment.imposedpressure import ImposedPressure
from xfv.src.cell.one_dimension_cell import OneDimensionCell
from xfv.src.data.data_container import DataContainer


class ImposedPressureTest(unittest.TestCase):
    """
    Test case utilisï¿½ pour test les fonctions du module 'Imposed Pressure'
    Marche bien mais appelle des mï¿½thodes de Cell et OneDimensionCell non vï¿½rifiï¿½es
    """

    def setUp(self) -> None:
        """
        Preparation of the test
        """
        # Creation of a fake DataContainer :
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_hydro.json")
        self.test_datacontainer = DataContainer(data_file_path)

        # Preparation of the test
        self._pressure = 0.
        self.my_imposed_pressure = ImposedPressure(self._pressure)
        self.cells = OneDimensionCell(1000)
        self.cells.pressure_field[:] = 1.
        self.ruptured_cells = np.ndarray([1000], dtype=np.bool, order='C')
        self.ruptured_cells[:] = False
        self.ruptured_cells[500] = True

    def tearDown(self) -> None:
        """
        End of the tests
        """
        DataContainer.clear()

    def test_apply_treatment(self) -> None:
        """
        Teste la mï¿½thode apply_treatment for ImposedPressure
        """
        self.my_imposed_pressure.apply_treatment(self.cells, self.ruptured_cells)
        self.cells.increment_variables()
        self.assertAlmostEqual(self.cells.pressure_field[500], self._pressure)


if __name__ == '__main__':
    unittest.main()
