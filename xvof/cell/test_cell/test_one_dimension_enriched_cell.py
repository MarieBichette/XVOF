#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module element1dupgraded
"""
import numpy as np
import unittest
import mock


from xvof.mesh.topology1d import Topology1D
from xvof.data.data_container import geometrical_props, material_props
# from xvof.data.data_container import numerical_props, properties
from xvof.data.data_container import numerical_props

from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.node.one_dimension_node import OneDimensionNode
from xvof.cell.one_dimension_enriched_cell import OneDimensionEnrichedCell


class OneDimensionEnrichedCellTest(unittest.TestCase):

    def setUp(self):
        """ Préparation des tests """
        num_props = numerical_props(0.2, 1.0, 0.35, 0.0)
        mat_props = material_props(1.0e+05, 0.0, 8129.,0., MieGruneisen(),"Enrichment", 0.5)
        geom_props = geometrical_props(1.0e-06)
        # props = properties(num_props, mat_props, geom_props)
        self.my_nodes = OneDimensionNode(2, poz_init=np.array([[-0.2,],[ 0.6,]]), vit_init=np.array([[0.,],[0.,]]), section=1.0e-06)
        #nodb = OneDimensionNode(1, poz_init=np.array([-0.2]), section=1.0e-06)
        #my_elem = el1d.OneDimensionCell(1)
         # attribué à my_elem:   props, 123, [noda, nodb])
        #self.my_elem_enr = el1denriched.OneDimensionEnrichedCell(my_elem, 0.5)
        self.my_elem_enr = OneDimensionEnrichedCell(1)
        self.topo1D = Topology1D(2, 1)



    def tearDown(self):
        pass


    def test_classical(self):
        """Test la propriété classical"""
        # Avant enrichissment :
        np.testing.assert_array_equal(self.my_elem_enr.classical, np.array([True]))
        # Après enrichissement : cell du milieu a été enrichie après application du critère
        self.my_elem_enr.classical[0] = False
        np.testing.assert_array_equal(self.my_elem_enr.classical, np.array([False]))
        # remise à zéros du paramètre
        self.my_elem_enr.classical[0] = True

    def test_enriched(self):
        """Test la propriété enriched"""
        # Avant enrichissement :
        np.testing.assert_array_equal(self.my_elem_enr.enriched, np.array([ False]))
        # Après enrichissement : cell du milieu a été enrichie après application du critère
        self.my_elem_enr.classical[0] = False
        np.testing.assert_array_equal(self.my_elem_enr.enriched, np.array([True]))
        # remise à zéros du paramètre
        self.my_elem_enr.classical[0] = True

    @mock.patch.object(Topology1D, 'nodes_belonging_to_cell', new_callable=mock.PropertyMock, return_value = np.array([0, 1]))
    @mock.patch.object(OneDimensionEnrichedCell, 'enriched', new_callable=mock.PropertyMock, return_value = np.array([True]))
    @mock.patch.object(OneDimensionEnrichedCell, 'left_size', new_callable=mock.PropertyMock, return_value = np.array([1.]))
    def test_get_left_part_coordinates(self, mock_left_size, mock_enriched, mock_connectivity):
        """Test de la methode get_left_part_coordinates"""
        nodes_coord = np.array([-0.2, 0.6])
        cell_left_coord = self.my_elem_enr.get_left_part_coordinates(self.topo1D, nodes_coord)
        np.testing.assert_array_equal(cell_left_coord, np.array([0.3]))

    @mock.patch.object(Topology1D, 'nodes_belonging_to_cell', new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    @mock.patch.object(OneDimensionEnrichedCell, 'enriched', new_callable=mock.PropertyMock, return_value=np.array([True]))
    @mock.patch.object(OneDimensionEnrichedCell, 'right_size', new_callable=mock.PropertyMock, return_value=np.array([1.]))
    def test_get_left_part_coordinates(self, mock_right_size, mock_enriched, mock_connectivity):
        """Test de la methode get_right_part_coordinates"""
        nodes_coord = np.array([-0.2, 0.6])
        cell_right_coord = self.my_elem_enr.get_left_part_coordinates(self.topo1D, nodes_coord)
        np.testing.assert_array_equal(cell_right_coord, np.array([0.1]))



    def test_compute_enriched_element_new_pressure(self):
        """Test de la methode compute_enriched_element_new_pressure"""
        # test sans external lib uniquement

        # un peu compliqué ...

    @mock.patch.object(OneDimensionEnrichedCell, "density", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionEnrichedCell, "left_size", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionEnrichedCell, "right_size", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionEnrichedCell, "enriched", new_callable= mock.PropertyMock,return_value= np.array([True]))
    def test_compute_enriched_element_new_density(self, mock_enriched, mock_left_size, mock_right_size, mock_density):
        """  Test of compute_enriched_element_new_density method    """
        self.my_elem_enr.density.current_left_value = np.array([1.])
        self.my_elem_enr.density.current_right_value = np.array([1.5])
        self.my_elem_enr.left_size.new_value = np.array([0.4])
        self.my_elem_enr.right_size.new_value = np.array([0.6])
        self.my_elem_enr.left_size.current_value = np.array([0.2])
        self.my_elem_enr.right_size.current_value = np.array([0.3])
        # import ipdb ; ipdb.set_trace()
        self.my_elem_enr.compute_enriched_elements_new_density()
        np.testing.assert_allclose(self.my_elem_enr.density.new_value, np.array([0.625]))
        np.testing.assert_allclose(self.my_elem_enr.density.new_enr_value, np.array([0.125]))



    @mock.patch.object(OneDimensionEnrichedCell, "density", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionEnrichedCell, "sound_velocity", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionEnrichedCell, "left_size", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionEnrichedCell, "right_size", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionEnrichedCell, "pseudo", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionEnrichedCell, "enriched", new_callable= mock.PropertyMock,return_value= np.array([True]))
    def test_compute_enriched_element_new_pseudo(self, mock_enriched, mock_pseudo, mock_right_size, mock_left_size,
                                                 mock_sound_velocity, mock_density):
        """Test de compute_enriched_element_new_pseudo method"""
        self.my_elem_enr.density.current_left_value = np.array([2000.])
        self.my_elem_enr.density.new_left_value = np.array([1000.])
        self.my_elem_enr.sound_velocity.current_left_value = np.array([0.])
        self.my_elem_enr.left_size.new_value = np.array([0.2])

        self.my_elem_enr.density.current_right_value = np.array([1000.])
        self.my_elem_enr.density.new_right_value = np.array([1000.])
        self.my_elem_enr.sound_velocity.current_right_value = np.array([0.])
        self.my_elem_enr.right_size.new_value = np.array([0.3])

        self.my_elem_enr.compute_enriched_elements_new_pseudo(1.)
        np.testing.assert_allclose(self.my_elem_enr.density.new_value, np.array([0.625]))
        np.testing.assert_allclose(self.my_elem_enr.density.new_enr_value, np.array([0.125])) #valeurs numériques à changer


    @mock.patch.object(Topology1D, 'nodes_belonging_to_cell', new_callable=mock.PropertyMock, return_value=np.array([0, 1]))
    def test_compute_enriched_elements_new_part_size(self, mock_connectivity):
        """ Test de la méthode OneDimensionEnrichedCell.compute_enriched_elements_new_part_size """
        vecteur_vitesse_noeuds = np.array([1.0, 1.0])
        vecteur_vitesse_enr_noeuds = np.array([-2.0, 2.0])
        self.my_elem_enr.compute_enriched_elements_new_part_size(1.0, mock_connectivity, vecteur_vitesse_enr_noeuds,
                                                                 vecteur_vitesse_noeuds)
        self.assertEqual(self.my_elem_enr.left_size.new_value, 0.4)
        self.assertEqual(self.my_elem_enr.right_size.new_value, 0.4)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
