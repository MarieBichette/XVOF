#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module node
"""
import xvof.node.node as nd
import numpy as np
import unittest


class NodeTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'Node'
    """
    def test_init(self):
        # L'index doit etre un entier
        with self.assertRaises(TypeError):
            nd.Node(index=3.0)
        # Vérification des valeurs par défaut
        default_node_2d = nd.Node(dim=2, index=1)
        np.testing.assert_array_equal(default_node_2d.coordt,
                                      np.zeros(2, dtype=float))
        np.testing.assert_array_equal(default_node_2d.coordtpdt,
                                      np.zeros(2, dtype=float))
        np.testing.assert_array_equal(default_node_2d.umundemi,
                                      np.zeros(2, dtype=float))
        np.testing.assert_array_equal(default_node_2d.upundemi,
                                      np.zeros(2, dtype=float))
        np.testing.assert_array_equal(default_node_2d.force,
                                      np.zeros(2, dtype=float))
        self.assertEqual(default_node_2d.masse, 0.)
        # Vérification de la cohérence entre la dimension du noeud et celle
        # des vecteurs positions et vitesses
        with self.assertRaises(SystemExit):
            node_1d_a = nd.Node(index=1, position_initiale=[0.1, 2.0])
        with self.assertRaises(SystemExit):
            node_1d_b = nd.Node(index=1, vitesse_initiale=[0.1, 2.0])
        # Vérification des valeurs passées en argument
        node_2d_a = nd.Node(dim=2, index=1, position_initiale=[0.1, -0.2],
                            vitesse_initiale=[-1.2, 2.0])
        np.testing.assert_array_equal(node_2d_a.coordt,
                                      np.array([0.1, -0.2]))
        np.testing.assert_array_equal(node_2d_a.coordtpdt,
                                      np.array([0.1, -0.2]))
        np.testing.assert_array_equal(node_2d_a.umundemi,
                                      np.array([-1.2, 2.0]))
        np.testing.assert_array_equal(node_2d_a.upundemi,
                                      np.array([-1.2, 2.0]))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
if __name__ == '__main__':
    unittest.main()