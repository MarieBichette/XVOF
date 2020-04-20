#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# pylint: disable=protected-access
"""
Classe de test du module node
"""
import unittest
import os
import numpy as np

import xfv.src.node.node as nd
from xfv.src.data.data_container import DataContainer


class NodeTest(unittest.TestCase):
    """
    Test case utilisï¿½ pour test les fonctions du module 'Node'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_hydro.json")
        self.test_datacontainer = DataContainer(data_file_path)

        self.vit_init = np.array([-1.5e+03, 1.2e+03, 0.3e+03], ndmin=2)
        self.poz_init = np.array([0.5, 0.025, -0.1], ndmin=2)
        # Crï¿½ation d'un noeud en 3D :
        self.my_node = nd.Node(1, position_initiale=self.poz_init,
                               vitesse_initiale=self.vit_init, dim=3)

    def tearDown(self):
        """
        Operations to be done after completing all the tests in the class
        """
        DataContainer.clear()

    def test_init(self):
        """
        Test du constructeur Node()
        """
        # Les valeurs par dï¿½faut de position et de vitesse doivent ï¿½tre des vecteurs nuls
        msg = "Par dï¿½faut, les vecteurs vitesses doivent ï¿½tre ï¿½gaux au vecteur nul!"
        msg2 = "Le vecteur force doit ï¿½tre initalisï¿½ au vecteur nul!"
        default_node_2d = nd.Node(1, np.array([0.5, 0.025], ndmin=2), dim=2)
        np.testing.assert_array_equal(default_node_2d.umundemi.flatten(),
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.upundemi.flatten(),
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.force, np.zeros([1,2], dtype=float), msg2)
        np.testing.assert_array_equal(default_node_2d.masse, np.array([[0., 0.], ]),
                                      "La masse doit ï¿½tre initalisï¿½e ï¿½ 0!")

        # Si la dimension du noeud ne correspond pas ï¿½ celle des vecteurs positions
        # et vitesses alors une exception de type SystemExit (pas top) doit ï¿½tre levï¿½e
        with self.assertRaises(SystemExit):
            nd.Node(1, position_initiale=np.array([[0.1, 2.0],]),
                    vitesse_initiale=np.array([0.1]), dim=1)
        with self.assertRaises(SystemExit):
            nd.Node(1, position_initiale=np.array([[0.1],]),
                    vitesse_initiale=np.array([0.1, 2.0]), dim=1)

        # Les vecteurs position et vitesse initiaux doivent ï¿½tre ceux passï¿½s en argument
        # du constructeur (sauf coordtpdt)
        node_2d_a = nd.Node(1, dim=2, position_initiale=np.array([[0.1, -0.2],]),
                            vitesse_initiale=np.array([[-1.2, 2.0],]))
        np.testing.assert_array_equal(node_2d_a.xt, np.array([[0.1, -0.2],]))
        np.testing.assert_array_equal(node_2d_a.xtpdt, np.array([[0.0, 0.0],]))
        np.testing.assert_array_equal(node_2d_a.umundemi, np.array([[-1.2, 2.0],]))
        np.testing.assert_array_equal(node_2d_a.upundemi, np.array([[-1.2, 2.0],]))

    def test_compute_new_coodinates(self):
        """
        Test de la mï¿½thode Node.compute_new_coodinates()
        """
        self.my_node.compute_new_coodinates(np.array([True]), delta_t=0.5e-06)
        np.testing.assert_allclose(self.my_node.xtpdt,
                                   np.array([0.49925, 0.0256, -0.09985], ndmin=2))

    def test_increment(self):
        """
        Test de la mï¿½thode Node.increment()
        """
        self.my_node.compute_new_coodinates(np.array([True]), delta_t=0.5e-06)
        self.my_node.increment()
        np.testing.assert_allclose(self.my_node.xt, np.array([0.49925, 0.0256, -0.09985], ndmin=2))
        np.testing.assert_array_equal(self.my_node.xt, self.my_node.xtpdt)
        np.testing.assert_array_equal(self.my_node.umundemi, self.my_node.upundemi)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
if __name__ == '__main__':
    unittest.main()
