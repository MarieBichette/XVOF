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
    def setUp(self):
        """
        Initialisation des tests
        """
        class Element():
            def __init__(self, masse):
                self.masse = masse

        ElemA = Element(3. / 4.)
        ElemB = Element(2. / 3.)
        vit_init = np.array([-1.5e+03, 1.2e+03, 0.3e+03])
        poz_init = np.array([0.5, 0.025, -0.1])
        self.my_node = nd.Node(dim=3, position_initiale=poz_init,
                               vitesse_initiale=vit_init)
        self.my_node.elements_voisins = [ElemA, ElemB]

    def test_init(self):
        """
        Test du constructeur Node()
        """
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
        msg = "Par défaut, les vecteurs positions et vitesses doivent "
        msg += "être égaux au vecteur nul!"
        msg2 = "Le vecteur force doit être initalisé au vecteur nul!"
        default_node_2d = nd.Node(dim=2, index=1)
        np.testing.assert_array_equal(default_node_2d.coordt,
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.coordtpdt,
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.umundemi,
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.upundemi,
                                      np.zeros(2, dtype=float), msg)
        np.testing.assert_array_equal(default_node_2d.force,
                                      np.zeros(2, dtype=float), msg2)
        self.assertEqual(default_node_2d.masse, 0., "La masse doit être \
        initalisée à 0!")
        #
        # Si la dimension du noeud ne correspond pas à celle
        # des vecteurs positions et vitesses alors une exception
        # de type SystemExit (pas top) doit être levée
        #
        with self.assertRaises(SystemExit):
            nd.Node(index=1, position_initiale=[0.1, 2.0])
        with self.assertRaises(SystemExit):
            nd.Node(index=1, vitesse_initiale=[0.1, 2.0])
        #
        # Les vecteurs position et vitesse initiaux doivent être
        # ceux passés en argument du constructeur (sauf coordtpdt)
        #
        node_2d_a = nd.Node(dim=2, index=1, position_initiale=[0.1, -0.2],
                            vitesse_initiale=[-1.2, 2.0])
        np.testing.assert_array_equal(node_2d_a.coordt,
                                      np.array([0.1, -0.2]))
        np.testing.assert_array_equal(node_2d_a.coordtpdt,
                                      np.array([0.0, 0.0]))
        np.testing.assert_array_equal(node_2d_a.umundemi,
                                      np.array([-1.2, 2.0]))
        np.testing.assert_array_equal(node_2d_a.upundemi,
                                      np.array([-1.2, 2.0]))

    def test_calculer_masse_wilkins(self):
        """
        Test de la fonction Node.calculer_masse_wilikins
        """
        self.my_node.calculer_masse_wilkins()
        self.assertEqual(self.my_node.masse, 0.7083333333333333)

    def test_calculer_nouvo_coord(self):
        """
        Test de la fonction Node.calculer_nouvo_coord
        """
        self.my_node.calculer_nouvo_coord(delta_t=0.5e-06)
        np.testing.assert_allclose(self.my_node.coordtpdt,
                                      np.array([0.49925, 0.0256, -0.09985]))

    def test_incrementer(self):
        """
        Test de la fonction Node.incrementer
        """
        self.my_node.calculer_nouvo_coord(delta_t=0.5e-06)
        self.my_node.incrementer()
        np.testing.assert_allclose(self.my_node.coordt,
                                      np.array([0.49925, 0.0256, -0.09985]))
        np.testing.assert_array_equal(self.my_node.coordt,
                                      self.my_node.coordtpdt)
        np.testing.assert_array_equal(self.my_node.umundemi,
                                      self.my_node.upundemi)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
if __name__ == '__main__':
    unittest.main()