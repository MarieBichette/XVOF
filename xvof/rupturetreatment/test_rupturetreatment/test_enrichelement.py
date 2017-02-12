#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module ImposedPressure
"""
import unittest
import xvof.rupturetreatment.imposedpressure as IP
import mock
import numpy as np

from xvof.cell.one_dimension_cell import OneDimensionCell


class ImposedPressureTest(unittest.TestCase):
    """
    Test case utilis� pour test les fonctions du module 'EnrichElement'
    """

    def setUp(self):
        '''
        Pr�paration du test
        '''


    def test_applyTreatment(self):
        '''
        Teste la m�thode applyTreatment for EnrichElement
        '''




if __name__ == '__main__':
    unittest.main()