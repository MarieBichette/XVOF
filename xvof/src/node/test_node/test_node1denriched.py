#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module OneDimensionEnrichedNode
"""
import numpy as np
import unittest
import mock

import xvof.src.node.one_dimension_enriched_node as nd1denr
import xvof.src.node.one_dimension_node as nd1d

from xvof.src.discontinuity.discontinuity import Discontinuity
from xvof.src.node.one_dimension_enriched_node import OneDimensionEnrichedNode
from xvof.src.node.one_dimension_node import OneDimensionNode
from xvof.src.mesh.topology1d import Topology1D


class OneDimensionEnrichedNodeTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'OneDimensionEnrichedNode'
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

        self.vit_init = np.zeros([6,1], dtype='float')
        self.vit_init[:,0] = [-1.5e+03, -0.6e+03, 1.2e+03, 0.0, 1.2e+03, 0.3e+03]
        self.poz_init = np.zeros([6,1], dtype='float')
        self.poz_init[:,0] = [0.6, 0.4, 0.35, 0.2, 0., -0.1]
        self.my_node_enr6 = nd1denr.OneDimensionEnrichedNode(6, self.poz_init, self.vit_init, section=1.0e-06)
        self.my_node_enr6._classical = np.array([True, True, True, True, True, True])

        self.vit_init = np.array([-1.5e+03], ndmin=2)
        self.poz_init = np.array([0.5], ndmin=2)
        self.my_node_enr = nd1denr.OneDimensionEnrichedNode(1, self.poz_init, self.vit_init, section=1.0e-06)
        self.my_node_enr.elements_voisins = [self.elem_droite, self.elem_gauche]

        discontinuity = Discontinuity(np.array([False, False, True, False, False, False]),
                                           np.array([False, False, False, True, False, False]))
        self.my_node_enr6.discontinuities = [discontinuity]




    def test_classical(self):
        """Test de la propriété classical du module OneDimensionEnrichedNode"""
        np.testing.assert_array_equal(self.my_node_enr6.classical, np.array([True, True, True, True, True,True]))
        # Après enrichissment : enrichissement des noeuds 1 et 2  (testé dans rupture_treatment.enrichment)
        self.my_node_enr6.classical[2] = False
        self.my_node_enr6.classical[3] = False
        np.testing.assert_array_equal(self.my_node_enr6.classical, np.array([True, True, False, False, True,True]))
        self.my_node_enr6.classical[2] = True
        self.my_node_enr6.classical[3] = True

    def test_enriched(self):
        """ Test de la propriété enriched du module OneDimensionEnrichedNode"""
        np.testing.assert_array_equal(self.my_node_enr6.enriched, np.array([False, False, False, False, False, False]))
        # Après enrichissment : enrichissement des noeuds 1 et 2  (testé dans rupture_treatment.enrichment)
        self.my_node_enr6.classical[2] = False
        self.my_node_enr6.classical[3] = False
        np.testing.assert_array_equal(self.my_node_enr6.enriched, np.array([False, False, True, True, False, False]))
        self.my_node_enr6.classical[2] = True
        self.my_node_enr6.classical[3] = True

    def test_enrichment_concerned(self):
        """ Test de la propriété enrichment_concerned du module OneDimensionEnrichedNode"""
        np.testing.assert_array_equal(self.my_node_enr6.enrichment_concerned, np.array([False, False, False, False, False, False]))
        # Après enrichissment : enrichissement des noeuds 1 et 2  (testé dans rupture_treatment.enrichment)
        self.my_node_enr6.classical[2] = False
        self.my_node_enr6.classical[3] = False
        np.testing.assert_array_equal(self.my_node_enr6.enrichment_concerned, np.array([False, True, True, True, True, False]))
        self.my_node_enr6.classical[2] = True
        self.my_node_enr6.classical[3] = True

    def test_enrichment_not_concerned(self):
        """ Test de la propriété enrichment_not_concerned du module OneDimensionEnrichedNode"""
        np.testing.assert_array_equal(self.my_node_enr6.enrichment_not_concerned, np.array([True, True, True, True, True, True]))
        # Après enrichissment : enrichissement des noeuds 1 et 2  (testé dans rupture_treatment.enrichment)
        self.my_node_enr6.classical[2] = False
        self.my_node_enr6.classical[3] = False
        np.testing.assert_array_equal(self.my_node_enr6.enrichment_not_concerned, np.array([True, False, False, False, False, True]))
        self.my_node_enr6.classical[2] = True
        self.my_node_enr6.classical[3] = True


    @mock.patch.object(OneDimensionEnrichedNode,'upundemi_enriched', new_callable=mock.PropertyMock,
                       return_value=np.array([0., 0., 1., 1., 0., 0.])) #verif signe
    def test_velocity_field(self, mock_up_enriched):
        self.my_node_enr6._upundemi = np.array([1., 1., 1., 1., 1., 1.])
        res_test = self.my_node_enr6.velocity_field
        # import ipdb; ipdb.set_trace()
        np.testing.assert_array_almost_equal(res_test, np.array([1., 1., 0., 2., 1., 1.]))


    def test_position_relative(self):
        """ Test de l'affectation de position_relative """
        #
        # Si position relative n'est pas -1 ou 1 une exception de type
        # SystemExit doit ëtre levée
        #
        with self.assertRaises(SystemExit):
            self.my_node_enr.position_relative = 2

    @mock.patch.object(Topology1D, "cells_in_contact_with_node", new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    @mock.patch.object(Topology1D, "getCellsInContactWithNode", new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    def test_enriched_nodes_compute_new_force(self, mock_get_cell, mock_cells_in_contact):
        """ Test de la méthode OneDimensionEnrichedNode.enriched_nodes_compute_new_force() """
        vecteur_pression_classique = np.array([self.elem_gauche.pression_t_plus_dt, self.elem_droite.pression_t_plus_dt])
        vecteur_pression_enrichie = np.array([self.elem_gauche._pression_t_plus_dt_enrichi, self.elem_droite._pression_t_plus_dt_enrichi])
        vecteur_pseudo_classique = np.array([self.elem_gauche.pseudo, self.elem_droite.pseudo])
        vecteur_pseudo_enrichie = np.array([self.elem_gauche._pseudo_plus_un_demi_enrichi, self.elem_droite._pseudo_plus_un_demi_enrichi])

        # Noeud à droite de la discontinuité
        self.my_node_enr.position_relative = -1
        self.my_node_enr.enriched_nodes_compute_new_force(Topology1D, vecteur_pression_classique, vecteur_pression_enrichie,
                                                            vecteur_pseudo_classique, vecteur_pseudo_enrichie)
        # force_classic testé dans one_dimension_node : np.testing.assert_array_equal(self.my_node.force_classique, np.array([2000.]))
        np.testing.assert_array_equal(self.my_node_enr._force_enriched, np.array([2750.]))
        # Noeud à gauche de la discontinuité
        self.my_node_enr.position_relative = 1
        self.my_node_enr.enriched_nodes_compute_new_force(Topology1D, vecteur_pression_classique, vecteur_pression_enrichie,
                                                            vecteur_pseudo_classique, vecteur_pseudo_enrichie)
        # force_classic testé dans one_dimension_node : np.testing.assert_array_equal(self.my_node.force_classique, np.array([2000.]))
        np.testing.assert_array_equal(self.my_node_enr._force_enriched, np.array([-3250.]))


    @mock.patch.object(OneDimensionEnrichedNode, 'upundemi_enriched',new_callable =mock.PropertyMock, return_value = np.array([1.]))
    @mock.patch.object(Discontinuity, 'mask_in_nodes', new_callable=mock.PropertyMock, return_value=np.array([True]))
    @mock.patch.object(Discontinuity, 'mask_out_nodes', new_callable=mock.PropertyMock, return_value=np.array([False]))
    def test_enriched_nodes_compute_new_coordinates_mask_in(self, mock_mask_out, mock_mask_in, mock_upenr):
        """Teste de la méthode enriched_nodes_compute_new_coordinates de la classe OneDimensionEnrichedNode"""
        self.my_node_enr._xtpdt = 0.
        self.my_node_enr.enriched_nodes_compute_new_coordinates(1.)
        np.testing.assert_equal(self.my_node_enr._xtpdt, -1.)

    @mock.patch.object(OneDimensionEnrichedNode, 'upundemi_enriched',new_callable =mock.PropertyMock, return_value = np.array([1.]))
    @mock.patch.object(Discontinuity, 'mask_in_nodes', new_callable=mock.PropertyMock, return_value=np.array([False]))
    @mock.patch.object(Discontinuity, 'mask_out_nodes', new_callable=mock.PropertyMock, return_value=np.array([True]))
    def test_enriched_nodes_compute_new_coordinates_mask_out(self, mock_mask_out, mock_mask_in, mock_upenr):
        """Teste de la méthode enriched_nodes_compute_new_coordinates de la classe OneDimensionEnrichedNode"""
        self.my_node_enr._xtpdt = 0.
        self.my_node_enr.enriched_nodes_compute_new_coordinates(1.)
        np.testing.assert_equal(self.my_node_enr._xtpdt, 1.)


    @mock.patch.object(Topology1D, "cells_in_contact_with_node", new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    @mock.patch.object(Topology1D, "getCellsInContactWithNode", new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    @mock.patch.object(OneDimensionEnrichedNode, 'umundemi_enriched',new_callable =mock.PropertyMock, return_value = np.array([1.]))
    def test_enriched_nodes_compute_new_velocity(self, mock_get_cell, mock_cells_in_contact):
        """ Test de la méthode OneDimensionEnrichedNode.enriched_nodes_compute_new_velocity() """
        vecteur_pression_classique = np.array([self.elem_gauche.pression_t_plus_dt, self.elem_droite.pression_t_plus_dt])
        vecteur_pression_enrichie = np.array([self.elem_gauche._pression_t_plus_dt_enrichi, self.elem_droite._pression_t_plus_dt_enrichi])
        vecteur_pseudo_classique = np.array([self.elem_gauche.pseudo, self.elem_droite.pseudo])
        vecteur_pseudo_enrichie = np.array([self.elem_gauche._pseudo_plus_un_demi_enrichi, self.elem_droite._pseudo_plus_un_demi_enrichi])

        #  Noeud à droite de la discontinuité
        self.my_node_enr.position_relative = -1
        self.my_node_enr.enriched_nodes_compute_new_force(Topology1D, vecteur_pression_classique, vecteur_pression_enrichie,
                                                            vecteur_pseudo_classique, vecteur_pseudo_enrichie)
        self.my_node_enr.enriched_nodes_compute_new_velocity(1.0e-01, np.array([True]),np.array([2.]))
        np.testing.assert_array_equal(self.my_node_enr.upundemi, np.array([400.0]))
        np.testing.assert_array_equal(self.my_node_enr.upundemi_enriched, np.array([550.0]))
        # Noeud à gauche de la discontinuité
        self.my_node_enr.position_relative = 1
        self.my_node_enr.enriched_nodes_compute_new_force(Topology1D, vecteur_pression_classique, vecteur_pression_enrichie,
                                                            vecteur_pseudo_classique, vecteur_pseudo_enrichie)
        self.my_node_enr.enriched_nodes_compute_new_velocity(1.0e-01, np.array([True]),np.array([2.]))
        np.testing.assert_array_equal(self.my_node_enr.upundemi, np.array([400.0]))
        np.testing.assert_array_equal(self.my_node_enr.upundemi_enriched, np.array([-650.0]))






if __name__ == '__main__':
    unittest.main()
