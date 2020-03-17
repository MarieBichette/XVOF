#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module node1d
"""
import numpy as np
import unittest
import mock

import xvof.src.node.one_dimension_node as nd1d
from xvof.src.mesh.topology1d import Topology1D
from xvof.src.node.one_dimension_node import OneDimensionNode


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
        self.my_node = nd1d.OneDimensionNode(1, np.array([0.],ndmin=2), np.array([0.],ndmin=2), section=1.0e-06)

    def test_elements_voisins(self):
        #  A laisser ici ? Pas de méthode elements_voisins dans node1d.
        #  Test <= 2 elements voisins implémentédanstopology1D
        #  Test à passer dans topologie1D

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
        self.my_node._elements_voisins = []
        self.my_node.elements_voisins = [self.elem_droite, self.elem_gauche]
        self.assertLessEqual(self.my_node.elements_voisins[0].coord,
                             self.my_node.elements_voisins[1].coord)

    @mock.patch.object(Topology1D, "cells_in_contact_with_node", new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    @mock.patch.object(Topology1D, "getCellsInContactWithNode", new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    def test_compute_new_force(self, mock_get_cell, mock_cells_in_contact):
        """ Test de la méthode Node1d.compute_new_force() """
        #
        # Test du calcul de la force en présence de 2 éléments voisins
        #
        pression_maille = np.array([self.elem_gauche.pression_t_plus_dt, self.elem_droite.pression_t_plus_dt])
        pseudo_maille = np.array([self.elem_gauche.pseudo, self.elem_droite.pseudo])
        self.my_node.elements_voisins = [self.elem_droite, self.elem_gauche]
        self.my_node.compute_new_force(Topology1D, pression_maille, pseudo_maille)
        np.testing.assert_array_equal(self.my_node.force, np.array([2000.]))
        self.my_node._elements_voisins = []
        #
        # Test du calcul de la force pour un noeud de bord droit
        #
        pression_maille = np.array([self.elem_gauche.pression_t_plus_dt])
        pseudo_maille = np.array([self.elem_gauche.pseudo])
        self.my_node.elements_voisins = [self.elem_gauche]
        self.my_node.compute_new_force(Topology1D, pression_maille, pseudo_maille)
        np.testing.assert_array_equal(self.my_node.force, np.array([3500.]))
        self.my_node._elements_voisins = []        
        #
        # Test du calcul de la force pour un noeud de bord gauche
        #
        pression_maille = np.array([self.elem_droite.pression_t_plus_dt])
        pseudo_maille = np.array([self.elem_droite.pseudo])
        self.my_node.elements_voisins = [self.elem_droite]
        self.my_node.compute_new_force(Topology1D, pression_maille, pseudo_maille)
        np.testing.assert_array_equal(self.my_node.force, np.array([-1500.]))

    @mock.patch.object(Topology1D, "cells_in_contact_with_node", new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    @mock.patch.object(Topology1D, "getCellsInContactWithNode", new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    def test_compute_new_velocity(self, mock_get_cell, mock_cells_in_contact):
        """ Test de la méthode Node1d.compute_new_velocity() """
        self.my_node.elements_voisins = [self.elem_droite, self.elem_gauche]
        # On a besoin de la force agissant sur le noeud
        pression_maille = np.array([self.elem_gauche.pression_t_plus_dt, self.elem_droite.pression_t_plus_dt])
        pseudo_maille = np.array([self.elem_gauche.pseudo, self.elem_droite.pseudo])
        self.my_node.compute_new_force(Topology1D, pression_maille, pseudo_maille)
        # On a besoin de la masse du noeud
        # Vérifiée par test dans test_node.py
        vecteur_masse_elements = np.array([self.elem_gauche.masse, self.elem_droite.masse])
        vecteur_nb_noeuds_par_element = np.array([1,1])
        masse = self.my_node.calculer_masse_wilkins(mock_get_cell, vecteur_masse_elements, vecteur_nb_noeuds_par_element)
        # Calcul de la vitesse (sans enrichisement)
        self.my_node.compute_new_velocity(1.0e-01, np.array([True]), masse)
        np.testing.assert_array_equal(self.my_node.upundemi, np.array([400.]))

    @mock.patch.object(OneDimensionNode, 'force', new_callable = mock.PropertyMock, return_value = np.array([-1500.]))
    @mock.patch.object(OneDimensionNode, 'section', new_callable = mock.PropertyMock, return_value = 10.)
    def test_apply_pressure(self, mock_section, mock_force):
        self.my_node.apply_pressure(0,100)
        np.testing.assert_array_equal(self.my_node.force, np.array([-500.]))
        #a revoir car fait appel à force dans le résultat=efface l'opération qui vient d'être faite

if __name__ == '__main__':
    unittest.main()
