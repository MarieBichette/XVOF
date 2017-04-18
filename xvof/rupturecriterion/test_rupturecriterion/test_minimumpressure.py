#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module node
"""
import numpy as np
import unittest

import mock

from xvof.cell.cell import Cell
from xvof.fields.field import Field
import xvof.rupturecriterion.minimumpressure as MinPressCrit

class MinimumPressureTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'ConstantPressure'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        self._prupturemin = 50.
        self._rupturecriterion = MinPressCrit.MinimumPressureCriterion(self._prupturemin)

        self.vecteur_pression_rupture = np.ones(1000) * 100
        self.vecteur_pression_rupture[500] = 25.

        self.vecteur_pression_pas_rupture = np.ones(1000) * 100
        self.vecteur_pression_pas_rupture[500] = 75.

        self.pressure_field = Field(1000, 0, 0)

        class CellTest() :
            self.pressure = self.pressure_field
            self.pressure.newvalue = self.vecteur_pression_pas_rupture

        self.Cell = CellTest()


    #
    # @mock.patch.object(Cell, "pressure", new_callable=mock.PropertyMock, return_value= pressure_field)
    # @mock.patch.object(Field, "new_value", new_callable=mock.PropertyMock, return_value = vecteur_pression_rupture)

    def test_checkCriterion_rupture(self, Cell):
        retour = self._rupturecriterion.checkCriterion(Cell(1000))
        self.assertFalse(retour[500])

# problème car fait appel à une autre méthode FieldManager() --> qu'est ce qu'on met dans le mock ???
        # 2 mock en argument sur la solution du haut....
#     @mock.patch.object(Cell, "pressure.new_value", new_callable=mock.PropertyMock, return_value=vecteur_pression_pas_rupture)
#     def test_checkCriterion_rupture(self, mock_connectivity):
#         retour = self._rupturecriterion.checkCriterion(Cell(1000))
#         self.assertTrue(retour[500])

if __name__ == '__main__':
    unittest.main()