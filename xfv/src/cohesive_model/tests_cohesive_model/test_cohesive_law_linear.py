#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Classe de test du module cohesive law
"""
import unittest
import os
import numpy as np

from xfv.src.cohesive_model.cohesive_law import CohesiveLaw
from xfv.src.data.data_container import DataContainer


class CohesiveLawLinearTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'CohesiveLaw' for linear cohesive law
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_hydro.xml")
        self.test_datacontainer = DataContainer(data_file_path)

    def tearDown(self):
        """
        Operations to be done after completing all the tests in the class
        """
        DataContainer.clear()




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
if __name__ == '__main__':
    unittest.main()
