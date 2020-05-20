# -*- coding: utf-8 -*-
# pylint: disable=protected-access, unused-argument
"""
Cell module unit tests
"""
import unittest
import unittest.mock as mock
import os
import numpy as np
from xfv.src.cell.one_dimension_enriched_cell_hansbo import OneDimensionHansboEnrichedCell
from xfv.src.cell.one_dimension_cell import OneDimensionCell
from xfv.src.data.data_container import DataContainer
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.utilities.stress_invariants_calculation import compute_second_invariant


class OneDimensionEnrichedHansboCellEPPTest(unittest.TestCase):
    """
    A class to test the OneDimensionHansboEnrichedCell module with enrichment,
    elasticity and plasticity
    """

    @classmethod
    def setUpClass(cls):
        """
        Tests setup for class
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_enrichment_epp.json")
        DataContainer(data_file_path)

    @classmethod
    def tearDownClass(cls):
        """
        Tear down after class tests
        """
        DataContainer.clear()
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

    def setUp(self):
        """
        Set up of tests
        """
        self.my_cells = OneDimensionHansboEnrichedCell(1)
        self.my_cells._classical = np.array([False])
        self.my_cells._external_library = None

        self.my_cells.energy.current_value = np.array([1.e+06])
        self.my_cells.pressure.current_value = np.array([1.5e+09])
        self.my_cells.density.current_value = np.array([8000.])
        self.my_cells.pseudo.current_value = np.array([1.e+08])
        self.my_cells.sound_velocity.current_value = np.array([300.])
        self.my_cells._deviatoric_stress_current = np.array([[5., 8., 3.]])

        self.my_cells.energy.new_value = np.array([0.8e+06])
        self.my_cells.pressure.new_value = np.array([1.3e+09])
        self.my_cells.density.new_value = np.array([8020.])
        self.my_cells.pseudo.new_value = np.array([1.e+08])
        self.my_cells.sound_velocity.new_value = np.array([302.])
        self.my_cells._deviatoric_stress_new = np.array([[4., 5., 6.]])
        self.my_cells._deviatoric_strain_rate = np.array([[1., 1., 1.]])

        self.my_cells.additional_dof_density.current_value = np.array([4000.])
        self.my_cells.additional_dof_density.new_value = np.array([4020.])
        self.my_cells.additional_dof_pressure.current_value = np.array([1.1e+09])
        self.my_cells.additional_dof_pressure.new_value = np.array([1.3e+09])
        self.my_cells.additional_dof_energy.current_value = np.array([1.e+06])
        self.my_cells.additional_dof_energy.new_value = np.array([0.8e+06])
        self.my_cells.additional_dof_artificial_viscosity.current_value = np.array([1.e+08])
        self.my_cells.additional_dof_artificial_viscosity.new_value = np.array([1.e+08])
        self.my_cells.additional_dof_sound_velocity.current_value = np.array([300.])
        self.my_cells.additional_dof_sound_velocity.new_value = np.array([302.])
        self.my_cells._additional_dof_deviatoric_stress_current = np.array([[3., 2., 1.],])
        self.my_cells._additional_dof_deviatoric_stress_new = np.array([[5., 12., 7.],])
        self.my_cells._additional_dof_deviatoric_stress_new = np.array([[5., 12., 7.], ])
        self.my_cells._additional_dof_deviatoric_strain_rate = np.array([[4., 3., 8.],])
        self.my_cells.additional_dof_yield_stress.current_value = np.array([10.])
        self.my_cells._additional_dof_equivalent_plastic_strain_rate = np.array([0.])
        self.my_cells._additional_dof_stress = np.array([[0., 0., 0.]])
        self.my_cells.left_part_size.current_value = np.array([0.2])
        self.my_cells.right_part_size.current_value = np.array([0.3])
        self.my_cells.left_part_size.new_value = np.array([0.4])
        self.my_cells.right_part_size.new_value = np.array([0.6])

        # configuration d'un mock 'discontinuity'
        config = {'mask_in_nodes': np.array([True, False]),
                  'mask_out_nodes': np.array([False, True]),
                  'position_in_ruptured_element':
                      DataContainer().material_target.failure_model.failure_treatment_value,
                  'mask_ruptured_cell': np.array([True]),
                  'ruptured_cell_id': 0,
                  'plastic_cells': np.array([False]),
                  }
        self.__patcher = mock.patch('xfv.src.discontinuity.discontinuity.Discontinuity',
                             spec=Discontinuity, **config)
        self.mock_disc = self.__patcher.start()

        self.test_data = DataContainer()  # pylint: disable=no-value-for-parameter

    def tearDown(self):
        self.__patcher.stop()
        return super().tearDown()

    def test_initialize_additional_dof(self):
        """
        Test of the method initialize_additional_dof
        """
        # todo : coder

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_enriched_pressure_field(self, mock_disc_list):
        """
        Test of the enriched field reconstruction
        (with the example of pressure field)
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        self.my_cells.additional_dof_pressure.current_value = np.array([1.e+09])
        self.my_cells.pressure._Field__current = np.array([1.])
        # le champ de pression est compos� de :
        # - valeur champ de pression � gauche pour cell 0
        # - valeur champ de pression � droite disc 0
        exact_field = np.array([1., 1.e+09])
        np.testing.assert_allclose(self.my_cells.pressure_field, exact_field)

    @mock.patch.object(OneDimensionCell, "apply_equation_of_state", spec=classmethod,
                       new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionCell, "add_elastic_energy_method", spec=classmethod,
                       new_callable=mock.MagicMock)
    def test_compute_enriched_elements_new_pressure_with_elasticity(self, mock_add_elasticity,
                                                                    mock_apply_eos):
        """
        Test of the method compute_enriched_elements_new_pressure pour Hansbo
        """
        # Configuration des mocks
        mock_apply_eos.return_value = self.my_cells.additional_dof_energy.new_value, \
                                      self.my_cells.additional_dof_pressure.new_value,\
                                      self.my_cells.additional_dof_sound_velocity.new_value
        mock_add_elasticity.return_value = self.my_cells.energy.new_value
        dt = 1  # pylint: disable=invalid-name

        self.my_cells.compute_enriched_elements_new_pressure(dt)

        mock_add_elasticity.assert_called()
        mock_add_elasticity.assert_called()
        mock_apply_eos.assert_any_call(
            self.my_cells, self.my_cells._target_eos,
            self.my_cells.additional_dof_density.current_value,
            self.my_cells.additional_dof_density.new_value,
            self.my_cells.additional_dof_pressure.current_value,
            self.my_cells.additional_dof_pressure.new_value,
            self.my_cells.additional_dof_energy.current_value,
            self.my_cells.additional_dof_energy.new_value,
            self.my_cells.additional_dof_artificial_viscosity.current_value,
            self.my_cells.additional_dof_sound_velocity.new_value)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionHansboEnrichedCell, "compute_discontinuity_borders_velocity",
                       spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_enriched_elements_new_part_size(self, mock_disc_border_velocity,
                                                     mock_disc_list):
        """
        Test of the method compute_enriched_elements_new_part_size pour
        OneDimensionEnrichedHansboCell
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        self.mock_disc.get_ruptured_cell_id.return_value = np.array([0])
        vecteur_vitesse_noeuds = np.array([[0.1, ], [0.15, ]])
        dt = 1.  # pylint: disable=invalid-name
        mock_disc_border_velocity.return_value = np.array([-0.05]), np.array([0.05])

        self.my_cells.compute_enriched_elements_new_part_size(dt, vecteur_vitesse_noeuds)
        np.testing.assert_allclose(self.my_cells.left_part_size.new_value, np.array([0.05]))
        np.testing.assert_allclose(self.my_cells.right_part_size.new_value, np.array([0.4]))

    def test_compute_enriched_element_new_density(self):
        """
        Test of the method compute_enriched_element_new_density method
        """
        self.my_cells.compute_enriched_elements_new_density()
        np.testing.assert_allclose(self.my_cells.density.new_value, np.array([4000.]))
        np.testing.assert_allclose(self.my_cells.additional_dof_density.new_value,
                                   np.array([2000.]))

    @mock.patch.object(OneDimensionCell, "compute_pseudo", spec=classmethod,
                       new_callable=mock.MagicMock)
    def test_compute_enriched_element_new_pseudo(self, mock_compute_pseudo):
        """
        Test of the method compute_enriched_element_new_pseudo method
        """
        mock_compute_pseudo.return_value = np.array([1.2])
        dt = 1.  # pylint: disable=invalid-name

        self.my_cells.compute_enriched_elements_new_pseudo(dt)

        mock_compute_pseudo.assert_any_call(
            dt, self.my_cells.density.current_value, self.my_cells.density.new_value,
            self.my_cells.left_part_size.new_value, self.my_cells.sound_velocity.current_value,
            self.test_data.numeric.a_pseudo, self.test_data.numeric.b_pseudo)

        mock_compute_pseudo.assert_any_call(
            dt, self.my_cells.additional_dof_density.current_value,
            self.my_cells.additional_dof_density.new_value,
            self.my_cells.right_part_size.new_value,
            self.my_cells.additional_dof_sound_velocity.current_value,
            self.test_data.numeric.a_pseudo, self.test_data.numeric.b_pseudo)

    @mock.patch.object(OneDimensionCell, "compute_time_step", spec=classmethod,
                       new_callable=mock.MagicMock)
    def test_compute_enriched_elements_new_time_step(self, mock_compute_dt):
        """
        Test of the method compute_enriched_elements_new_time_step
        """
        mock_compute_dt.return_value = np.array([0.])

        self.my_cells.compute_enriched_elements_new_time_step()

        mock_compute_dt.assert_any_call(self.test_data.numeric.cfl,
                                        self.test_data.numeric.cfl_pseudo,
                                        self.my_cells.density.current_value,
                                        self.my_cells.density.new_value,
                                        self.my_cells.left_part_size.new_value,
                                        self.my_cells.sound_velocity.new_value,
                                        self.my_cells.pseudo.current_value,
                                        self.my_cells.pseudo.new_value)
        mock_compute_dt.assert_called_with(
            self.test_data.numeric.cfl, self.test_data.numeric.cfl_pseudo,
            self.my_cells.additional_dof_density.current_value,
            self.my_cells.additional_dof_density.new_value,
            self.my_cells.right_part_size.new_value,
            self.my_cells.additional_dof_sound_velocity.new_value,
            self.my_cells.additional_dof_artificial_viscosity.current_value,
            self.my_cells.additional_dof_artificial_viscosity.new_value)

    def test_compute_enriched_stress_tensor(self):
        """
        Test of the method compute_enriched_stress_tensor
        """
        exact_stress_droite = np.array([[1., -1., 2.]])
        # Left part : computed in the classical method

        # Test for the left part of enriched cells
        self.my_cells._additional_dof_deviatoric_stress_new = np.array([[2., 0., 3.]])
        self.my_cells.additional_dof_pressure.new_value = np.array([1.])
        self.my_cells.additional_dof_artificial_viscosity.new_value = np.array([0.])
        self.my_cells.compute_enriched_stress_tensor()
        np.testing.assert_allclose(self.my_cells.additional_dof_stress, exact_stress_droite)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionHansboEnrichedCell, "_compute_discontinuity_borders_velocity",
                       spec=classmethod, new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionCell, "general_method_deviator_strain_rate",
                       spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_enriched_deviatoric_strain_rate(self, mock_compute_d, mock_disc_borders,
                                                     mock_disc_list):
        """
        Test of the method compute_enriched_deviatoric_strain_rate
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        dt = 1.  # pylint: disable=invalid-name
        node_coord_new = np.array([[0.,], [1.,]])
        node_velocity_new = np.array([[-1, ], [1., ]])
        Discontinuity.additional_dof_velocity_new = np.array([[[1., ], [3., ]]])
        u1d_arr = Discontinuity.additional_dof_velocity_new[:, 1]
        u2g_arr = Discontinuity.additional_dof_velocity_new[:, 0]
        Discontinuity.discontinuity_position = np.array([0.5])
        Discontinuity.in_nodes = np.array([0])
        Discontinuity.out_nodes = np.array([1])
        u_disc_g = np.array([-0.5])
        u_disc_d = np.array([0.5])
        mock_disc_borders.return_value = u_disc_g, u_disc_d
        mock_compute_d.return_value = np.array([[1., 1., 1.], [1., 1., 1.]])
        self.my_cells.plastic_enr_cells = np.array([True])

        self.my_cells.compute_enriched_deviatoric_strain_rate(dt, node_coord_new, node_velocity_new)

        mock_disc_borders.assert_called_with(self.mock_disc, node_velocity_new)
        mock_compute_d.assert_called()

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionCell, "general_method_deviator_strain_rate", spec=classmethod,
                       new_callable=mock.MagicMock)
    def test_compute_enriched_deviatoric_stress_tensor(self,
                                                       mock_compute_d, mock_disc_list):
        """
        Test of the method compute_enriched_deviatoric_stress_tensor
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        self.my_cells._classical = np.array([False])
        dt = 1.  # pylint: disable=invalid-name
        coord_noeud_new = np.array([[-1.], [0, ]])
        vitesse_noeud_new = np.array([[50, ], [-20, ]])
        mock_compute_d.return_value = np.array([[2., -1, -1], [2., -1, -1]])
        Discontinuity.additional_dof_velocity_new = np.array([[[1., ], [3., ]]])

        self.my_cells.additional_dof_shear_modulus.new_value = np.array([14.])
        self.my_cells._additional_dof_deviatoric_stress_current = np.array([[0., 0., 0.]])
        self.my_cells._additional_dof_deviatoric_strain_rate = np.array([[2., -1, -1]])
        exact_new_s_right = np.array([[56., -28, -28]])

        self.my_cells.compute_enriched_deviatoric_stress_tensor(coord_noeud_new,
                                                                vitesse_noeud_new, dt)
        mock_compute_d.assert_called()
        np.testing.assert_allclose(self.my_cells._additional_dof_deviatoric_stress_new,
                                   exact_new_s_right)

    def test_compute_enriched_shear_modulus(self):
        """
        Test of the method compute_enriched_shear_modulus
        """
        target_model = self.test_data.material_target
        self.my_cells.compute_enriched_shear_modulus(
            target_model.constitutive_model.elasticity_model.build_shear_modulus_obj())

        np.testing.assert_allclose(self.my_cells.additional_dof_shear_modulus.new_value,
                                   target_model.initial_values.shear_modulus_init)

    def test_compute_enriched_equivalent_plastic_strain_rate(self):
        """
        Test of the method compute_enriched_equivalent_plastic_strain_rate
        """
        dt = 1.  # pylint: disable=invalid-name
        exact_equivalent_plastic_strain_rate_right = np.array([0.33333333333333331])
        deviatoric_stress_new = np.array([[4, -2., -2.]])
        j2 = np.sqrt(compute_second_invariant(deviatoric_stress_new))
        shear_modulus = np.array([4.])
        yield_stress = np.array([2.])

        res = self.my_cells._compute_equivalent_plastic_strain_rate(j2, shear_modulus, yield_stress, dt)

        np.testing.assert_allclose(res, exact_equivalent_plastic_strain_rate_right)

    def test_apply_plastic_correction_on_enriched_deviatoric_stress_tensor(self):
        """
        Test of the method apply_plastic_correction_on_enriched_deviatoric_stress_tensor
        """
        deviatoric_stress_new = np.array([[5., -2.5, -2.5], ])
        j2 = np.sqrt(compute_second_invariant(deviatoric_stress_new))
        yield_stress = np.array([10.])
        exact_deviator_plastic_right = np.array([[6.666667, -3.333333, -3.333333]])
        radial_return = self.my_cells._compute_radial_return(j2, yield_stress)
        deviatoric_stress_new *= radial_return[np.newaxis].T
        np.testing.assert_allclose(deviatoric_stress_new, exact_deviator_plastic_right, rtol=1.e-5)

    def test_compute_enriched_yield_stress(self):
        """
        Test of the method compute_enriched_shear_modulus
        """
        target_model = self.test_data.material_target
        self.my_cells.compute_enriched_yield_stress(
            target_model.constitutive_model.plasticity_model.build_yield_stress_obj())

        np.testing.assert_allclose(self.my_cells.additional_dof_yield_stress.new_value,
                                   target_model.initial_values.yield_stress_init)


if __name__ == "__main__":
    unittest.main()
