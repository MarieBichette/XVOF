#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module Topology1d
"""
import unittest
import os
from xfv.src.mesh.topology1d import Topology1D
from xfv.src.data.data_container import DataContainer


class Topology1dTest(unittest.TestCase):
    """
    Test case utilis� pour test les fonctions du module 'Topology1d'
    """
    def setUp(self):
        """
        Preparation des tests
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_hydro.json")
        self.test_datacontainer = DataContainer(data_file_path)
        self.topology = Topology1D(4, 3)

    def tearDown(self):
        DataContainer.clear()

    def test_elements_voisins(self):
        """
        Test de Node1d.elements_voisins
        """
        #
        # En 1D affecter plus de deux elements � un noeud doit lever
        # une exception de type SystemExit
        #
        with self.assertRaises(IndexError):
            self.topology.add_cell_in_contact_with_node(1, 2)


if __name__ == '__main__':
    unittest.main()
