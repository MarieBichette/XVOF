#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Classe de test du module cohesive zone model
"""
import unittest
import os
import numpy as np

from xfv.src.cohesive_model.linear_cohesive_law import LinearCohesiveZoneModel
from xfv.src.data.data_container import DataContainer


class LinearCohesiveZoneModelTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'LinearCohesiveZoneModel'
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
