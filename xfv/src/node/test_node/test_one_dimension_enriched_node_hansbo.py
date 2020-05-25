# -*- coding: utf-8 -*-
# pylint: disable=protected-access, unused-argument
"""
Classe de test du module OneDimensionEnrichedNode
"""
import unittest
import unittest.mock as mock
import os
import numpy as np

from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.node.one_dimension_enriched_node_hansbo import OneDimensionHansboEnrichedNode
from xfv.src.data.data_container import DataContainer


class OneDimensionEnrichedNodeHansboTest(unittest.TestCase):
    """
    Test case used for module 'OneDimensionHansboEnrichedNode'
    """
    def setUp(self):
        """
        Preparation of the unit tests
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_hydro.json")
        self.test_datacontainer = DataContainer(data_file_path)

        self.vit_init = np.zeros([2, 1], dtype='float')
        self.vit_init[:, 0] = [1.2e+03, 0.0]
        self.poz_init = np.zeros([2, 1], dtype='float')
        self.poz_init[:, 0] = [1., 2.]
        self.my_nodes = OneDimensionHansboEnrichedNode(2, self.poz_init, self.vit_init,
                                                       section=1.0e-06)
        self.my_nodes._classical = np.array([False, False])

        # configuration d'un mock 'discontinuity'
        config = {'mask_in_nodes': np.array([True, False]),
                  'mask_out_nodes': np.array([False, True]),
                  'label': 1,
                  'position_in_ruptured_element':
                      self.test_datacontainer.material_target.failure_model.failure_treatment_value,
                  'mask_ruptured_cell': np.array([True]),
                  'ruptured_cell_id': np.array([0]),
                  'plastic_cells': np.array([False]),
                  'additional_dof_force': np.array([[1., ], [2., ]]),
                  '_additional_dof_velocity_new': np.zeros([2, 1])
                  }
        patcher = mock.patch('xfv.src.discontinuity.discontinuity.Discontinuity',
                             spec=Discontinuity, **config)
        self.mock_discontinuity = patcher.start()

    def tearDown(self):
        """
        Operations to be done after completing all the tests in the class
        """
        DataContainer.clear()

    def test_classical(self):
        """
        Test de la propri�t� classical du module OneDimensionEnrichedNode
        """
        np.testing.assert_array_equal(self.my_nodes.classical, np.array([False, False]))

    def test_enriched(self):
        """
        Test de la propri�t� enriched du module OneDimensionEnrichedNode
        """
        np.testing.assert_array_equal(self.my_nodes.enriched, np.array([True, True]))

    def test_enrichment_concerned(self):
        """
        Test de la propri�t� enrichment_concerned du module
        """
        np.testing.assert_array_equal(self.my_nodes.enrichment_concerned, np.array([True, True]))

    def test_enrichment_not_concerned(self):
        """
        Test de la propri�t� enrichment_not_concerned du module
        """
        np.testing.assert_array_equal(self.my_nodes.enrichment_not_concerned,
                                      np.array([False, False]))

    def test_compute_complete_velocity_field(self):
        """
        Test de la m�thode compute_complete_velocity_field de la classe
        OneDimensionHansboEnrichedNodes
        """
        self.my_nodes._upundemi = np.array([1., 1.])
        self.my_nodes.compute_complete_velocity_field()
        np.testing.assert_array_almost_equal(self.my_nodes.velocity_field, np.array([1., 1.]))

    def test_compute_additional_dof_new_velocity(self):
        """
        Test of the method compute_additional_dof_new_velocity
        """
        # Discontinuity.discontinuity_list.return_value = [self.mock_discontinuity]
        # inv_mass_additional = np.array([2., 2.])
        # self.mock_discontinuity.__class__.additional_dof_velocity_new = np.array([[[0., ], [0., ]]])
        # self.mock_discontinuity.__class__.additional_dof_velocity_current = np.array([[[1., ], [3., ]]])
        # self.mock_discontinuity.__class__.additional_dof_force = np.array([[[1., ], [1., ]]])
        # self.my_nodes.compute_additional_dof_new_velocity(1., inv_mass_additional)
        # np.testing.assert_array_almost_equal(self.mock_discontinuity.__class__.additional_dof_velocity_new[0],
        self.mock_discontinuity.additional_dof_velocity_new = np.array([[0., ], [0., ]])
        self.mock_discontinuity.mass_matrix_enriched = mock.PropertyMock()
        self.mock_discontinuity.mass_matrix_enriched.inverse_enriched_mass_matrix_enriched_dof = np.array([[2., ], [2., ]])
        self.mock_discontinuity.additional_dof_velocity_current = np.array([[1., ], [3., ]])
        self.mock_discontinuity.additional_dof_force = np.array([[1., ], [1., ]])
        self.my_nodes.compute_additional_dof_new_velocity(self.mock_discontinuity, 1.)
        np.testing.assert_array_almost_equal(self.mock_discontinuity.additional_dof_velocity_new,
                                             np.array([[3., ], [5., ]]))

    def test_coupled_enrichment_terms_compute_new_velocity(self):
        """
        Test of the method coupled_terms_compute_new_velocity
        """
        self.mock_discontinuity.additional_dof_velocity_new = np.array([[0., ], [0., ]])
        self.mock_discontinuity.mass_matrix_enriched = mock.PropertyMock()
        self.mock_discontinuity.mass_matrix_enriched.inverse_enriched_mass_matrix_coupling_dof = np.array([[1., 2.], [2., 1.]])
        self.my_nodes._force = np.array([[1., ], [1., ]])
        self.my_nodes._upundemi = np.array([[1., ], [1., ]])
        self.my_nodes.coupled_enrichment_terms_compute_new_velocity(self.mock_discontinuity, 1.)

        np.testing.assert_array_equal(self.my_nodes._upundemi, np.array([[6., ], [5., ]]))
        np.testing.assert_array_equal(self.mock_discontinuity.additional_dof_velocity_new,
                                      np.array([[3., ], [3., ]]))

    def test_enriched_nodes_compute_new_coordinates(self):
        """
        Test de la méthode enriched_nodes_compute_new_coordinates de la classe
        """
        self.mock_discontinuity.additional_dof_coordinates_current = np.array([[1., ], [3., ]])
        self.mock_discontinuity.additional_dof_coordinates_new = np.array([[-1., ], [-3., ]])
        self.mock_discontinuity.additional_dof_velocity_new = np.array([[-3., ], [4., ]])
        self.my_nodes.enriched_nodes_compute_new_coordinates(self.mock_discontinuity, 1.)
        np.testing.assert_array_equal(self.mock_discontinuity._additional_dof_coordinates_new,
                                      np.array([[-2., ], [7., ]]))

    def test_reinitialize_kinematics_after_contact(self):
        """
        Test of the method reinitialize_kinematics_after_contact
        """
        self.mock_discontinuity.mask_disc_nodes = np.ones([2], dtype=bool)
        self.my_nodes._umundemi = np.array([[-0.5, ], [1.5, ]])
        self.my_nodes._xt = np.array([[0., ], [1., ]])
        self.my_nodes._upundemi = np.array([[-1.5, ], [2.5, ]])
        self.my_nodes._xtpdt = np.array([[2., ], [3., ]])
        self.my_nodes.reinitialize_kinematics_after_contact(self.mock_discontinuity)
        np.testing.assert_array_equal(self.my_nodes._upundemi, np.array([[-0.5, ], [1.5, ]]))
        np.testing.assert_array_equal(self.my_nodes._xtpdt, np.array([[0., ], [1., ]]))

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_enriched_nodes_new_force(self, mock_disc_list):
        """
        Test de la m�thode enriched_nodes_compute_new_force
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_discontinuity]
        self.mock_discontinuity.position_in_ruptured_element = 0.25
        self.mock_discontinuity.get_ruptured_cell_id = 0
        contrainte_classique = np.array([2.])
        contrainte_enr = np.array([2.])
        self.my_nodes._force = np.array([[4., ], [2., ]])
        self.my_nodes._section = 1.
        self.my_nodes.compute_enriched_nodes_new_force(contrainte_classique, contrainte_enr)

        np.testing.assert_array_almost_equal(
            self.mock_discontinuity.additional_dof_force, np.array([[-0.5, ], [1.5, ]]))
        np.testing.assert_almost_equal(self.my_nodes._force, np.array([[5.5, ], [1.5, ]]))

    def test_apply_force_on_discontinuity_boundaries(self):
        """
        Test of the method apply_force_on_discontinuity_boundaries
        """
        stress = 50.
        self.my_nodes._section = 1.
        self.my_nodes._force = np.array([[200., ], [400., ]])
        self.mock_discontinuity.position_in_ruptured_element = 0.5
        self.mock_discontinuity.additional_dof_force = np.array([[-100., ], [300., ]])
        self.my_nodes.apply_force_on_discontinuity_boundaries(self.mock_discontinuity, stress)
        np.testing.assert_almost_equal(self.my_nodes._force, np.array([[225., ], [375., ]]))
        np.testing.assert_array_almost_equal(
            self.mock_discontinuity.additional_dof_force, np.array([[-75., ], [275., ]]))

    @unittest.skip("Mod�le coh�sif pas revu")
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_nodes_cohesive_forces(self, mock_disc_list):
        """
        Test de la m�thode compute_enriched_nodes_cohesive_forces
        """
        # Test des autres cas : la discontinuit� est en train de s'ouvrir
        Discontinuity.discontinuity_list.return_value = [self.mock_discontinuity]
        self.my_nodes._force = np.array([[0., ], [0., ]])
        self.mock_discontinuity.position_in_ruptured_element = 0.25
        self.mock_discontinuity.additional_dof_force = np.array([[0., ], [0., ]])
        self.mock_discontinuity.discontinuity_opening = 0.5
        self.mock_discontinuity.cohesive_force.new_value = [0.]

        exact_force_classic = np.array([[7.5, ], [-2.5, ]])
        exact_force_enriched = np.array([[2.5, ], [-7.5, ]])

        # d�finition de ouverture old
        self.my_nodes._xt = np.array([[0., ], [1., ]])
        self.mock_discontinuity.right_part_size.current_value = np.array([0.4])
        self.mock_discontinuity.left_part_size.current_value = np.array([0.4])
        xd_old = self.my_nodes.xt[self.mock_discontinuity.mask_out_nodes] \
                 - self.mock_discontinuity.right_part_size.current_value
        xg_old = self.my_nodes.xt[self.mock_discontinuity.mask_in_nodes] + \
                 self.mock_discontinuity.left_part_size.current_value
        ouverture_ecaille_old = (xd_old - xg_old)[0][0]
        np.testing.assert_allclose(ouverture_ecaille_old, np.array([0.2]))

        # definition de ouverture new
        self.my_nodes._xtpdt = np.array([[0., ], [1., ]])
        self.mock_discontinuity.right_part_size.new_value = np.array([0.3])
        self.mock_discontinuity.left_part_size.new_value = np.array([0.3])
        xd_old = self.my_nodes.xtpdt[self.mock_discontinuity.mask_out_nodes] - \
                 self.mock_discontinuity.right_part_size.new_value
        xg_old = self.my_nodes.xtpdt[self.mock_discontinuity.mask_in_nodes] + \
                 self.mock_discontinuity.left_part_size.new_value
        ouverture_ecaille_new = (xd_old - xg_old)[0][0]
        np.testing.assert_allclose(ouverture_ecaille_new, np.array([0.4]))

        self.my_nodes.compute_enriched_nodes_cohesive_forces()

        np.testing.assert_allclose(self.my_nodes.force, exact_force_classic)
        np.testing.assert_allclose(self.mock_discontinuity.additional_dof_force,
                                   exact_force_enriched)

if __name__ == '__main__':
    unittest.main()
