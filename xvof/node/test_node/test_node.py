#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module node
"""
import numpy as np
import unittest
import mock

import xvof.node.node as nd
from xvof.mesh.topology1d import Topology1D


class NodeTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'Node'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        class Element():
            def __init__(self, masse):
                self.masse = masse

        self.ElemA = Element(3. / 4.)
        self.ElemB = Element(2. / 3.)
        self.vit_init = np.array([-1.5e+03, 1.2e+03, 0.3e+03], ndmin=2)
        self.poz_init = np.array([0.5, 0.025, -0.1], ndmin=2)
        self.my_node = nd.Node(1, position_initiale=self.poz_init,
                               vitesse_initiale=self.vit_init, dim=3) # crée 2 noeuds Node en dimension 3
        self.my_node.elements_voisins = [self.ElemA, self.ElemB]

    def test_init(self):
        """ Test du constructeur Node() """
        #
        # Si l'index n'est pas un entier alors
        # une exception de type TypeError doit être
        # levée
        #
        with self.assertRaises(TypeError):
            nd.Node(index=3.0)
        #
        # Les valeurs par défaut de position et de vitesse
        # doivent être des vecteurs nuls
        #
        msg = "Par défaut, les vecteurs positions et vitesses doivent être égaux au vecteur nul!"
        msg2 = "Le vecteur force doit être initalisé au vecteur nul!"
        default_node_2d = nd.Node(dim=2, index=1)
        np.testing.assert_array_equal(default_node_2d.xt,
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.xtpdt,
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.umundemi,
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.upundemi,
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.force,
                                      np.zeros(2, dtype=float), msg2)
        self.assertEqual(default_node_2d.masse, 0., "La masse doit être initalisée à 0!")
        #
        # Si la dimension du noeud ne correspond pas à celle
        # des vecteurs positions et vitesses alors une exception
        # de type SystemExit (pas top) doit être levée
        #
        with self.assertRaises(SystemExit):
            nd.Node(1, position_initiale=[0.1, 2.0], vitesse_initiale=[0.1], dim=1)
        with self.assertRaises(SystemExit):
            nd.Node(1, position_initiale=[0.1], vitesse_initiale=[0.1, 2.0],dim=1)
        #
        # Les vecteurs position et vitesse initiaux doivent être
        # ceux passés en argument du constructeur (sauf coordtpdt)
        #
        node_2d_a = nd.Node(1, dim=2, position_initiale=[0.1, -0.2],
                            vitesse_initiale=[-1.2, 2.0])
        np.testing.assert_array_equal(node_2d_a.xt,
                                      np.array([0.1, -0.2]))
        np.testing.assert_array_equal(node_2d_a.xtpdt,
                                      np.array([0.0, 0.0]))
        np.testing.assert_array_equal(node_2d_a.umundemi,
                                      np.array([-1.2, 2.0]))
        np.testing.assert_array_equal(node_2d_a.upundemi,
                                      np.array([-1.2, 2.0]))


    @mock.patch.object(Topology1D, 'getCellsInContactWithNode', new_callable = mock.PropertyMock, return_value = np.array([0, 1]))
    def test_calculer_masse_wilkins(self, mock_topologie):
        """ Test de la méthode Node.calculer_masse_wilikins() """
        vecteur_masse_elements = np.array([self.ElemA.masse, self.ElemB.masse])
        vecteur_nb_noeuds_par_element = np.array([1,1])
        self.my_node.calculer_masse_wilkins(mock_topologie, vecteur_masse_elements, vecteur_nb_noeuds_par_element)
        self.assertEqual(self.my_node.masse, 0.7083333333333333)

    def test_compute_new_coodinates(self):
        """ Test de la méthode Node.compute_new_coodinates() """
        self.my_node.compute_new_coodinates(delta_t=0.5e-06)
        np.testing.assert_allclose(self.my_node.xtpdt,
                                      np.array([0.49925, 0.0256, -0.09985], ndmin=2))

    def test_increment(self):
        """ Test de la méthode Node.increment() """
        self.my_node.compute_new_coodinates(delta_t=0.5e-06)
        self.my_node.increment()
        np.testing.assert_allclose(self.my_node.xt,
                                      np.array([0.49925, 0.0256, -0.09985], ndmin=2))
        np.testing.assert_array_equal(self.my_node.xt,
                                      self.my_node.xtpdt)
        np.testing.assert_array_equal(self.my_node.umundemi,
                                      self.my_node.upundemi)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
if __name__ == '__main__':
    unittest.main()
