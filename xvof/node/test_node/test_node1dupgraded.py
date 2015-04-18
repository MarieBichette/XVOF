#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module node1dupgraded
"""
import unittest

import numpy as np
import xvof.node.node1d as nd1d
import xvof.node.node1dupgraded as nd1dup


class Node1dUpgradedTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'Node1dUpgraded'
    """
    def setUp(self):
        """
        Préparation des tests unitaires
        """
        class element():
            def __init__(self, poz, pressure, pseudo, masse):
                self.coord = poz
                self.pression_t_plus_dt = pressure
                self._pression_t_plus_dt_enrichi = pressure * 0.5
                self.pseudo = pseudo
                self._pseudo_plus_un_demi_enrichi = pseudo * 0.5
                self.masse = masse
        self.elem_gauche = element(np.array([-0.5]), 2.5e+09, 1.0e+09, 3. / 4.)
        self.elem_droite = element(np.array([0.5]), 1.0e+09, 0.5e+09, 1. / 4.)
        my_node = nd1d.Node1d(1, section=1.0e-06)
        my_node.elements_voisins = [self.elem_droite, self.elem_gauche]
        self.my_node = nd1dup.Node1dUpgraded(my_node)
        self.my_node.calculer_masse_wilkins()

    def test_position_relative(self):
        """ Test de l'affectation de position_relative """
        #
        # Si position relative n'est pas -1 ou 1 une exception de type
        # SystemExit doit ëtre levée
        #
        with self.assertRaises(SystemExit):
            self.my_node.position_relative = 2

    def test_calculer_nouvo_force(self):
        """ Test de la méthode Node1dUpgraded.calculer_nouvo_force() """
        # Noeud à droite de la discontinuité
        self.my_node.position_relative = -1
        self.my_node.calculer_nouvo_force()
        np.testing.assert_array_equal(self.my_node.force_classique,
                                      np.array([2000.]))
        np.testing.assert_array_equal(self.my_node.force_enrichi,
                                      np.array([2750.]))
        # Noeud à gauche de la discontinuité
        self.my_node.position_relative = 1
        self.my_node.calculer_nouvo_force()
        np.testing.assert_array_equal(self.my_node.force_classique,
                                      np.array([2000.]))
        np.testing.assert_array_equal(self.my_node.force_enrichi,
                                      np.array([-3250.]))

    def test_calculer_nouvo_vitesse(self):
        """ Test de la méthode Node1dUpgraded.calculer_nouvo_vitesse() """
        #  Noeud à droite de la discontinuité
        self.my_node.position_relative = -1
        self.my_node.calculer_nouvo_force()
        self.my_node.calculer_nouvo_vitesse(1.0e-01)
        np.testing.assert_array_equal(self.my_node.upundemi_classique,
                                      np.array([400.0]))
        np.testing.assert_array_equal(self.my_node.upundemi_enrichi,
                                      np.array([550.0]))
        # Noeud à gauche de la discontinuité
        self.my_node.position_relative = 1
        self.my_node.calculer_nouvo_force()
        self.my_node.calculer_nouvo_vitesse(1.0e-01)
        np.testing.assert_array_equal(self.my_node.upundemi_classique,
                                      np.array([400.0]))
        np.testing.assert_array_equal(self.my_node.upundemi_enrichi,
                                      np.array([-650.0]))

if __name__ == '__main__':
    unittest.main()
