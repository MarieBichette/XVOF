#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module one_dimension_enriched_cell_Moes
"""
import numpy as np
import unittest
import mock
from xvof.mesh.topology1d import Topology1D
from xvof.cell.one_dimension_enriched_cell_Moes import OneDimensionMoesEnrichedCell
from xvof.data.data_container import DataContainer
from xvof.discontinuity.discontinuity import Discontinuity


class OneDimensionEnrichedMoesCellTest(unittest.TestCase):

    def setUp(self):
        """ Préparation des tests """
        data_file_path = "//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml"
        self.test_datacontainer = DataContainer(data_file_path)
        self.my_elements = OneDimensionMoesEnrichedCell(1)

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_element_new_pressure(self, mock_disc_list):
        """Test de la methode compute_enriched_element_new_pressure
        Without external lib """

        self.my_elements._classical = np.array([False])
        self.my_elements.energy.current_value = np.array([1.e+06])
        self.my_elements.pressure.current_value = np.array([1.5e+09])

        self.my_elements.density.current_value = np.array([8000.])
        self.my_elements.pseudo.current_value = np.array([1.e+08])
        self.my_elements.sound_velocity.current_value = np.array([300.])

        self.my_elements.energy.new_value = np.array([0.8e+06])
        self.my_elements.pressure.new_value = np.array([1.3e+09])
        self.my_elements.density.new_value = np.array([8020.])
        self.my_elements.pseudo.new_value = np.array([1.e+08])
        self.my_elements.sound_velocity.new_value = np.array([302.])

        disc = mock.MagicMock(Discontinuity)
        type(disc).mask_in_nodes = mock.PropertyMock(return_value=np.array([True, False]))
        type(disc).mask_out_nodes = mock.PropertyMock(return_value=np.array([False, True]))
        type(disc).position_in_ruptured_element = mock.PropertyMock(return_value=0.25)
        Discontinuity.discontinuity_list.return_value = [disc]

        disc.additional_dof_density = mock.PropertyMock()
        disc.additional_dof_density.current_value = np.array([8000.])
        disc.additional_dof_density.new_value = np.array([8100.])

        disc.additional_dof_pressure = mock.PropertyMock()
        disc.additional_dof_pressure.current_value = np.array([2.e+09])
        disc.additional_dof_pressure.new_value = np.array([1.e+09])

        disc.additional_dof_energy = mock.PropertyMock()
        disc.additional_dof_energy.current_value = np.array([8.])
        disc.additional_dof_energy.new_value = np.array([8.2])

        disc.additional_dof_artificial_viscosity = mock.PropertyMock()
        disc.additional_dof_artificial_viscosity.current_value = np.array([1.e+09])

        disc.additional_dof_sound_velocity = mock.PropertyMock()
        disc.additional_dof_sound_velocity.current_value = np.array([300.])

        self.my_elements._external_library = None
        self.my_elements.compute_enriched_elements_new_pressure()

        np.testing.assert_array_almost_equal(self.my_elements.pressure.new_value, np.array([0.]))
        np.testing.assert_array_almost_equal(self.my_elements.energy.new_value, np.array([0.]))
        np.testing.assert_array_almost_equal(self.my_elements.sound_velocity.new_value, np.array([0.]))
        np.testing.assert_array_almost_equal(disc.additional_dof_pressure.new_value, np.array([0.]))
        np.testing.assert_array_almost_equal(disc.additional_dof_energy.new_value, np.array([0.]))
        np.testing.assert_array_almost_equal(disc.additional_dof_sound_velocity.new_value, np.array([0.]))

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_elements_new_part_size(self, mock_disc_list):
        """ Test de la méthode compute_enriched_elements_new_part_size pour OneDimensionEnrichedMoesCell"""
        disc = mock.MagicMock(Discontinuity)
        type(disc).mask_in_nodes = mock.PropertyMock(return_value=np.array([True, False]))
        type(disc).mask_out_nodes = mock.PropertyMock(return_value=np.array([False, True]))
        type(disc).position_in_ruptured_element = mock.PropertyMock(return_value=0.25)
        type(disc).additional_dof_velocity_new = mock.PropertyMock(return_value=np.array([[1., ], [3., ]]))
        Discontinuity.discontinuity_list.return_value = [disc]

        disc.left_part_size = mock.PropertyMock()
        disc.left_part_size.current_value = np.array([1.0])

        disc.right_part_size = mock.PropertyMock()
        disc.right_part_size.current_value = np.array([1.0])

        vecteur_vitesse_noeuds = np.array([[1., ], [1.5, ]])
        self.my_elements.compute_enriched_elements_new_part_size(1., vecteur_vitesse_noeuds)
        # problèmes dimensions ???!!!
        np.testing.assert_array_almost_equal(disc.left_part_size.new_value, np.array([[0.625]]))
        np.testing.assert_array_almost_equal(disc.right_part_size.new_value, np.array([[2.875]]))

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_element_new_density(self, mock_disc_list):
        """
        Test of compute_enriched_element_new_density method
        """
        topo = mock.MagicMock(Topology1D)
        type(topo).cells_in_contact_with_node = mock.PropertyMock(return_value=np.array([[-1, 0], [0, 1]]))

        disc = mock.MagicMock(Discontinuity)
        type(disc).mask_in_nodes = mock.PropertyMock(return_value=np.array([True, False]))
        type(disc).mask_out_nodes = mock.PropertyMock(return_value=np.array([False, True]))
        type(disc).position_in_ruptured_element = mock.PropertyMock(return_value=0.25)
        Discontinuity.discontinuity_list.return_value = [disc]

        disc.left_part_size = mock.PropertyMock()
        disc.left_part_size.current_value = np.array([0.2])
        disc.left_part_size.new_value = np.array([0.4])

        disc.right_part_size = mock.PropertyMock()
        disc.right_part_size.current_value = np.array([0.3])
        disc.right_part_size.new_value = np.array([0.6])

        disc.additional_dof_density = mock.PropertyMock()
        disc.additional_dof_density.current_value = np.array([1.5])

        self.my_elements.density.current_value = np.array([1.])

        self.my_elements.compute_enriched_elements_new_density(topo)
        np.testing.assert_allclose(self.my_elements.density.new_value, np.array([0.5]))
        np.testing.assert_allclose(disc.additional_dof_density.new_value, np.array([0.75]))

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_element_new_pseudo(self, mock_disc_list):
        """
        Test de compute_enriched_element_new_pseudo method
        """
        topo = mock.MagicMock(Topology1D)
        type(topo).cells_in_contact_with_node = mock.PropertyMock(return_value=np.array([[-1, 0], [0, 1]]))

        self.my_elements.density.new_value = np.array([3.5])
        self.my_elements.density.current_value = np.array([2.])
        self.my_elements.sound_velocity.current_value = np.array([1.])

        disc = mock.MagicMock(Discontinuity)
        type(disc).mask_in_nodes = mock.PropertyMock(return_value=np.array([True, False]))
        Discontinuity.discontinuity_list.return_value = [disc]
        disc.additional_dof_density = mock.PropertyMock()
        disc.additional_dof_density.current_value = np.array([0.])
        disc.additional_dof_density.new_value = np.array([-0.5])
        disc.right_part_size = mock.PropertyMock()
        disc.right_part_size.new_value = np.array([0.6])
        disc.additional_dof_sound_velocity = mock.PropertyMock()
        disc.additional_dof_sound_velocity.current_value = np.array([0.])

        disc.left_part_size = mock.PropertyMock()
        disc.left_part_size.new_value = np.array([0.4])

        type(DataContainer().numeric).a_pseudo = mock.PropertyMock(return_value=1.5)
        type(DataContainer().numeric).b_pseudo = mock.PropertyMock(return_value=0.2)

        self.my_elements.compute_enriched_elements_new_pseudo(1., topo)
        np.testing.assert_allclose(self.my_elements.pseudo.new_value, np.array([0.37461333333]))
        np.testing.assert_allclose(disc.additional_dof_artificial_viscosity.new_value, np.array([-0.0520533333]))

    @unittest.skipUnless(DataContainer("//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml").
                         material_target.failure_model.type_of_enrichment == "Moes", "Moes method not maintained")
    def test_compute_enriched_elements_new_time_step(self):
        """
        Test de la méthode compute_enriched_elements_new_time_step
        """
        # à coder quand on aura corrigé la méthode de calcul

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()