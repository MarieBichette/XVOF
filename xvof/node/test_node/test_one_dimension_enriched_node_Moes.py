#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module OneDimensionEnrichedNode
"""
import numpy as np
import unittest
import mock
from xvof.discontinuity.discontinuity import Discontinuity
from xvof.node.one_dimension_enriched_node_Moes import OneDimensionMoesEnrichedNode
from xvof.mesh.topology1d import Topology1D
from xvof.data.data_container import DataContainer


class OneDimensionEnrichedNodeMoesTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'OneDimensionEnrichedMoesNode'
    """
    def setUp(self):
        """
        Préparation des tests unitaires
        """
        class element:
            def __init__(self, poz, pressure, pseudo, masse):
                self.coord = poz
                self.pression_new = pressure
                self.pseudo = pseudo
                self.masse = masse
        self.elem_0 = element(np.array([-0.5]), 2.5e+09, 1.0e+09, 1. / 4.)
        self.elem_1 = element(np.array([0.5]), 1.0e+09, 0.5e+09, 1. / 4.)
        self.elem_2 = element(np.array([1.5]), 2.0e+09, 0.e+09, 1. / 2.)

        self.vit_init = np.zeros([4, 1], dtype='float')
        self.vit_init[:, 0] = [-1.5e+03, 1.2e+03, 0.0, 0.3e+03]
        self.poz_init = np.zeros([4, 1], dtype='float')
        self.poz_init[:, 0] = [0., 1., 2., 3.]
        self.my_nodes = OneDimensionMoesEnrichedNode(4, self.poz_init, self.vit_init, section=1.0e-06)
        self.my_nodes._classical = np.array([True, False, False, True])

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    def test_classical(self):
        """
        Test de la propriété classical du module OneDimensionEnrichedMoesNode
        """
        np.testing.assert_array_equal(self.my_nodes.classical, np.array([True, False, False, True]))
        self.my_nodes._classical[1] = True
        self.my_nodes._classical[2] = True
        np.testing.assert_array_equal(self.my_nodes.classical, np.array([True, True, True, True]))
        self.my_nodes._classical[1] = False
        self.my_nodes._classical[2] = False

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    def test_enriched(self):
        """
        Test de la propriété enriched du module OneDimensionEnrichedMoesNode
        """
        np.testing.assert_array_equal(self.my_nodes.enriched, np.array([False, True, True, False]))
        self.my_nodes._classical[1] = True
        self.my_nodes._classical[2] = True
        np.testing.assert_array_equal(self.my_nodes.enriched, np.array([False, False, False, False]))
        self.my_nodes._classical[1] = False
        self.my_nodes._classical[2] = False

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    def test_enrichment_concerned(self):
        """
        Test de la propriété enrichment_concerned du module OneDimensionEnrichedMoesNode
        """
        np.testing.assert_array_equal(self.my_nodes.enrichment_concerned, np.array([True, True, True, True]))
        self.my_nodes._classical[1] = True
        self.my_nodes._classical[2] = True
        np.testing.assert_array_equal(self.my_nodes.enrichment_concerned, np.array([False, False, False, False]))
        self.my_nodes._classical[1] = False
        self.my_nodes._classical[2] = False

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    def test_enrichment_not_concerned(self):
        """
        Test de la propriété enrichment_not_concerned du module OneDimensionEnrichedMoesNode
        """
        np.testing.assert_array_equal(self.my_nodes.enrichment_not_concerned, np.array([False, False, False, False]))
        self.my_nodes.classical[1] = True
        self.my_nodes.classical[2] = True
        np.testing.assert_array_equal(self.my_nodes.enrichment_not_concerned, np.array([True, True, True, True]))
        self.my_nodes.classical[1] = False
        self.my_nodes.classical[2] = False

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_complete_velocity_field(self, mock_disc_list):
        """
        Test de la méthode compute_complete_velocity_field de la classe OneDimensionMoesEnrichedNodes
        """
        disc = mock.MagicMock(Discontinuity)
        type(disc).mask_in_nodes= mock.PropertyMock(return_value=np.array([False, True, False, False]))
        type(disc).mask_out_nodes = mock.PropertyMock(return_value=np.array([False, False, True, False]))
        type(disc)._additional_dof_velocity_new = mock.PropertyMock(return_value=np.array([1., 1.]))
        Discontinuity.discontinuity_list.return_value = [disc]

        self.my_nodes._upundemi = np.array([1., 1., 1., 1.])
        self.my_nodes.compute_complete_velocity_field()
        np.testing.assert_array_almost_equal(self.my_nodes.velocity_field, np.array([1., 0., 2., 1.]))

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_coupled_enrichment_terms_compute_new_velocity(self, mock_disc_list):
        disc = mock.MagicMock(Discontinuity)
        type(disc).mask_in_nodes = mock.PropertyMock(return_value=np.array([False, True, False, False]))
        type(disc).mask_out_nodes = mock.PropertyMock(return_value=np.array([False, False, True, False]))
        type(disc)._additional_dof_velocity_new = mock.PropertyMock(return_value=np.zeros([2, 1]))
        type(disc).additional_dof_force = mock.PropertyMock(return_value=np.array([[1., ], [2., ]]))
        Discontinuity.discontinuity_list.return_value = [disc]

        inv_masse_couplage = np.array([[1., 2.], [2., 1.], [0., 3.], [3., 0]])
        self.my_nodes._force = np.array([[1., ], [1., ], [1., ], [1., ]])
        self.my_nodes._upundemi = np.array([[1., ], [1., ], [1., ], [1., ]])
        self.my_nodes.coupled_enrichment_terms_compute_new_velocity(1., inv_masse_couplage)

        np.testing.assert_array_equal(self.my_nodes._upundemi,np.array([[6., ], [5., ],  [7., ], [4., ]]))
        np.testing.assert_array_equal(disc._additional_dof_velocity_new, np.array([[6., ], [6., ]]))

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_enriched_nodes_compute_new_coordinates(self, mock_disc_list):
        disc = mock.MagicMock(Discontinuity)
        type(disc).mask_in_nodes = mock.PropertyMock(return_value=np.array([False, True, False, False]))
        type(disc).mask_out_nodes = mock.PropertyMock(return_value=np.array([False, False, True, False]))
        type(disc).additional_dof_velocity_new = mock.PropertyMock(return_value=np.ones([2, 1]))
        Discontinuity.discontinuity_list.return_value = [disc]

        self.my_nodes._xtpdt = np.array([[0., ], [1., ], [2., ], [3., ]])
        self.my_nodes.enriched_nodes_compute_new_coordinates(0.5)

        np.testing.assert_array_equal(self.my_nodes.xtpdt, np.array([[0., ], [0.5, ], [2.5, ], [3., ]]))

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_enriched_nodes_compute_new_force(self, mock_disc_list):
        # Mock disc :
        disc = mock.MagicMock(Discontinuity)
        type(disc).mask_in_nodes = mock.PropertyMock(return_value=np.array([False, True, False, False]))
        type(disc).mask_out_nodes = mock.PropertyMock(return_value=np.array([False, False, True, False]))
        type(disc).additional_dof_force = mock.PropertyMock(return_value=np.array([[1., ], [2., ]]))
        type(disc.additional_dof_pressure).new_value = mock.PropertyMock(return_value=np.array([2.e+09]))
        type(disc.additional_dof_artificial_viscosity).new_value = mock.PropertyMock(return_value=np.array([1.e+09]))
        type(disc).position_in_ruptured_element = mock.PropertyMock(return_value=0.1)
        Discontinuity.discontinuity_list.return_value = [disc]
        # Mock topo :
        topo = mock.MagicMock(Topology1D)
        type(topo).cells_in_contact_with_node = mock.PropertyMock(
            return_value=np.array([[-1, 0], [0, 1], [1, 2], [2, -1]]))
        topo.getCellsInContactWithNode.side_effect = [np.array([0, -1]), np.array([2, -1])]

        vecteur_pression_classique = np.array(
            [self.elem_0.pression_new, self.elem_1.pression_new, self.elem_2.pression_new])
        vecteur_pseudo_classique = np.array([self.elem_0.pseudo, self.elem_1.pseudo, self.elem_2.pseudo])
        self.my_nodes._force = np.array([[1., ], [1., ], [1., ], [1., ]])
        self.my_nodes._section = 1.e-09

        self.my_nodes.enriched_nodes_compute_new_force(topo, vecteur_pression_classique, vecteur_pseudo_classique)

        np.testing.assert_almost_equal(disc.additional_dof_force, np.array([[-7.7, ], [2.2, ]]))
        np.testing.assert_almost_equal(self.my_nodes._force, np.array([[1.], [-1.4], [3.4], [1.]]))


if __name__ == '__main__':
    unittest.main()