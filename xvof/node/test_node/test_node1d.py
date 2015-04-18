#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module node1d
"""
import unittest

import numpy as np
import xvof.node.node1d as nd1d


class Node1dTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'Node1d'
    """
    def setUp(self):
        """
        Préparation des tests
        """
        class element():
            def __init__(self, poz, pressure, pseudo, masse):
                self.coord = poz
                self.pression_t_plus_dt = pressure
                self.pseudo = pseudo
                self.masse = masse
        self.elem_gauche = element(np.array([-0.5]), 2.5e+09, 1.0e+09, 3. / 4.)
        self.elem_droite = element(np.array([0.5]), 1.0e+09, 0.5e+09, 1. / 4.)
        self.elem_nul = element(np.array([0.15]), 3.3e+09, 0.0, 0.0)
        self.my_node = nd1d.Node1d(1, section=1.0e-06)

    def test_elements_voisins(self):
        """ Test de Node1d.elements_voisins = """
        #
        # En 1D affecter plus de deux éléments à un noeud doit lever
        # une exception de type SystemExit
        #
        with self.assertRaises(SystemExit):
            self.my_node.elements_voisins = [self.elem_droite,
                                             self.elem_gauche, self.elem_nul]
        #
        # Les éléments doivent être triés selon les coordonnées croissantes
        #
        self.my_node.elements_voisins = [self.elem_droite, self.elem_gauche]
        self.assertLessEqual(self.my_node.elements_voisins[0].coord,
                             self.my_node.elements_voisins[1].coord)

    def test_calculer_nouvo_force(self):
        """ Test de la méthode Node1d.calculer_nouvo_force() """
        #
        # Test du calcul de la force en présence de 2 éléments voisins
        #
        self.my_node.elements_voisins = [self.elem_droite, self.elem_gauche]
        self.my_node.calculer_nouvo_force()
        np.testing.assert_array_equal(self.my_node.force, np.array([2000.]))
        #
        # Test du calcul de la force pour un noeud de bord droit
        #
        self.my_node.elements_voisins = [self.elem_gauche]
        self.my_node.calculer_nouvo_force()
        np.testing.assert_array_equal(self.my_node.force, np.array([3500.]))
        #
        # Test du calcul de la force pour un noeud de bord gauche
        #
        self.my_node.elements_voisins = [self.elem_droite]
        self.my_node.calculer_nouvo_force()
        np.testing.assert_array_equal(self.my_node.force, np.array([-1500.]))

    def test_calculer_nouvo_vitesse(self):
        """ Test de la méthode Node1d.calculer_nouvo_vitesse() """
        self.my_node.elements_voisins = [self.elem_droite, self.elem_gauche]
        # On a besoin de la force agissant sur le noeud
        self.my_node.calculer_nouvo_force()
        # On a besoin de la masse du noeud
        # Vérifiée par test dans test_node.py
        self.my_node.calculer_masse_wilkins()
        # Calcul de la vitesse
        self.my_node.calculer_nouvo_vitesse(1.0e-01)
        np.testing.assert_array_equal(self.my_node.upundemi, np.array([400.]))

if __name__ == '__main__':
    unittest.main()
