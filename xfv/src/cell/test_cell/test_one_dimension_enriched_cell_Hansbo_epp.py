#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Cell module unit tests
"""
import numpy as np
import unittest
import mock
import os
from xfv.src.mesh.topology1d import Topology1D
from xfv.src.cell.one_dimension_enriched_cell_Hansbo import OneDimensionHansboEnrichedCell
from xfv.src.cell.one_dimension_enriched_cell import OneDimensionEnrichedCell
from xfv.src.cell.one_dimension_cell import OneDimensionCell
from xfv.src.data.data_container import DataContainer
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.rheology.constantshearmodulus import ConstantShearModulus


class OneDimensionEnrichedHansboCellEPPTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Tests setup for class
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_enrichment_epp.xml")
        DataContainer(data_file_path)

    @classmethod
    def tearDownClass(cls):
        DataContainer.clear()
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
        pass

    def setUp(self):
        """
        Pr�paration des tests
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

        # configuration d'un mock 'discontinuity'
        config = {'mask_in_nodes': np.array([True, False]),
                  'mask_out_nodes': np.array([False, True]),
                  'position_in_ruptured_element':
                      DataContainer().material_target.failure_model.failure_treatment_value,
                  'mask_ruptured_cell': np.array([True]),
                  'ruptured_cell_id': np.array([0]),
                  'plastic_cells': False,
                  'additional_dof_density.current_value':np.array([4000.]),
                  'additional_dof_density.new_value': np.array([4020.]),
                  'additional_dof_pressure.current_value': np.array([1.1e+09]),
                  'additional_dof_pressure.new_value': np.array([1.3e+09]),
                  'additional_dof_energy.current_value': np.array([1.e+06]),
                  'additional_dof_energy.new_value': np.array([0.8e+06]),
                  'additional_dof_artificial_viscosity.current_value': np.array([1.e+08]),
                  'additional_dof_artificial_viscosity.new_value': np.array([1.e+08]),
                  'additional_dof_sound_velocity.current_value': np.array([300.]),
                  'additional_dof_sound_velocity.new_value': np.array([302.]),
                  'additional_dof_deviatoric_stress_current': np.array([[3., 2., 1.],]),
                  'additional_dof_deviatoric_stress_new': np.array([[5., 12., 7.],]),
                  '_additional_dof_deviatoric_stress_new': np.array([[5., 12., 7.], ]),
                  'additional_dof_deviatoric_strain_rate': np.array([[4., 3., 8.],]),
                  'additional_dof_yield_stress.current_value': np.array([10.]),
                  '_additional_dof_equivalent_plastic_strain_rate': np.array([0.]),
                  'additional_dof_velocity_new': np.array([[1., ], [3., ]]),
                  'additional_dof_stress': np.array([[0., 0., 0.]]),
                  'left_part_size.current_value': np.array([0.2]),
                  'right_part_size.current_value': np.array([0.3]),
                  'left_part_size.new_value': np.array([0.4]),
                  'right_part_size.new_value': np.array([0.6]),
                  }
        patcher = mock.patch('xfv.src.discontinuity.discontinuity.Discontinuity',
                             spec=Discontinuity, **config)
        self.mock_disc = patcher.start()

    def tearDown(self):
        pass

    def test_initialize_additional_dof(self):
        """ Test la m�thode initialize_additional_dof"""
        # todo : coder



    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_enriched_pressure_field(self, mock_disc_list):
        """
        Test de la reconstruction des champs enrichis
        On prend l'exemple de la pression mais le test est valable pour tous les champs.
        """
        # configuration de la premi�re discontinuit�
        config = {'label': 1,
                  'ruptured_cell_id': np.array([0]),
                  'additional_dof_pressure.current_value': np.array([1.e+09])}
        patcher = mock.patch('xfv.src.discontinuity.discontinuity.Discontinuity',
                             spec=Discontinuity, **config)
        discontinuity_1 = patcher.start()
        # configuration de la deuxi�me discontinuit�
        config = {'label': 2,
                  'ruptured_cell_id': np.array([1]),
                  'additional_dof_pressure.current_value': np.array([5.e+09])}
        patcher = mock.patch('xfv.src.discontinuity.discontinuity.Discontinuity',
                             spec=Discontinuity, **config)
        discontinuity_2 = patcher.start()

        # on a besoin de faire comme si on avait 2 cells classiques qui rompent
        self.my_cells.pressure._Field__current = np.array([1., 2.])

        Discontinuity.discontinuity_list.return_value = [discontinuity_2, discontinuity_1]
        # le champ de pression est compos� de :
        # - valeur champ de pression � gauche pour cell 0
        # - valeur champ de pression � droite disc 0
        # - valeur champ de pression � gauche pour cell 1
        # - valeur champ de pression "� droite" pour disc 1
        exact_field = np.array([1., 1.e+09, 2., 5.e+09])
        np.testing.assert_allclose(self.my_cells.pressure_field, exact_field)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionCell, "apply_equation_of_state", spec=classmethod,
                       new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionCell, "add_elastic_energy_method", spec=classmethod,
                       new_callable=mock.MagicMock)
    def test_compute_enriched_elements_new_pressure_with_elasticity(self, mock_add_elasticity,
                                                                    mock_apply_eos, mock_disc_list):
        """
        Test de la m�thode compute_enriched_elements_new_pressure pour Hansbo
        """
        # Configuration des mocks
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        mock_apply_eos.return_value = self.my_cells.energy.new_value, \
                                      self.my_cells.pressure.new_value,\
                                      self.my_cells.sound_velocity.new_value
        mock_add_elasticity.return_value = self.my_cells.energy.new_value
        dt = 1

        self.my_cells.compute_enriched_elements_new_pressure(dt)

        mock_add_elasticity.assert_called()
        mock_add_elasticity.assert_called()
        mock_apply_eos.assert_any_call(self.my_cells,
                                       DataContainer().material_target.constitutive_model.eos,
                                       self.my_cells.density.current_value,
                                       self.my_cells.density.new_value,
                                       self.my_cells.pressure.current_value,
                                       self.my_cells.pressure.new_value,
                                       self.my_cells.energy.current_value,
                                       self.my_cells.energy.new_value,
                                       self.my_cells.pseudo.current_value,
                                       self.my_cells.sound_velocity.new_value)
        mock_apply_eos.assert_any_call(self.my_cells,
                                       DataContainer().material_target.constitutive_model.eos,
                                       self.mock_disc.additional_dof_density.current_value,
                                       self.mock_disc.additional_dof_density.new_value,
                                       self.mock_disc.additional_dof_pressure.current_value,
                                       self.mock_disc.additional_dof_pressure.new_value,
                                       self.mock_disc.additional_dof_energy.current_value,
                                       self.mock_disc.additional_dof_energy.new_value,
                                       self.mock_disc.additional_dof_artificial_viscosity.current_value,
                                       self.mock_disc.additional_dof_sound_velocity.new_value)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionHansboEnrichedCell, "compute_discontinuity_borders_velocity",
                       spec=classmethod, new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionEnrichedCell, "compute_new_left_right_size",
                       spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_enriched_elements_new_part_size(self, mock_compute_size,
                                                     mock_disc_border_velocity,
                                                     mock_disc_list):
        """
        Test de la m�thode compute_enriched_elements_new_part_size pour
        OneDimensionEnrichedHansboCell
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        vecteur_vitesse_noeuds = np.array([[1., ], [1.5, ]])
        dt = 1.
        u1g = vecteur_vitesse_noeuds[self.mock_disc.mask_in_nodes]
        u2d = vecteur_vitesse_noeuds[self.mock_disc.mask_out_nodes]
        ug = mock.PropertyMock(return_value=np.array([5.]))
        ud = mock.PropertyMock(return_value=np.array([6.]))
        mock_disc_border_velocity.return_value = ug, ud

        self.my_cells.compute_enriched_elements_new_part_size(dt, vecteur_vitesse_noeuds)
        mock_disc_border_velocity.assert_called_with(self.mock_disc,
                                                     vecteur_vitesse_noeuds)
        mock_compute_size.assert_called_with(dt, self.mock_disc, u1g, u2d, ug, ud)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionEnrichedCell, "compute_new_left_right_density",
                       spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_enriched_element_new_density(self, mock_compute_density, mock_disc_list):
        """
        Test of compute_enriched_element_new_density method
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        mock_compute_density.return_value = np.array([2.]), np.array([3.])

        self.my_cells.compute_enriched_elements_new_density()
        mock_compute_density.assert_called_with(self.my_cells.density.current_value,
                                                self.mock_disc.additional_dof_density.current_value,
                                                self.mock_disc)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionCell, "compute_pseudo", spec=classmethod,
                       new_callable=mock.MagicMock)
    def test_compute_enriched_element_new_pseudo(self, mock_compute_pseudo, mock_disc_list):
        """
        Test de compute_enriched_element_new_pseudo method
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        mock_compute_pseudo.return_value = np.array([1.2])
        dt = 1.

        self.my_cells.compute_enriched_elements_new_pseudo(dt)

        mock_compute_pseudo.assert_any_call(dt, self.my_cells.density.current_value,
                                            self.my_cells.density.new_value,
                                            self.mock_disc.left_part_size.new_value,
                                            self.my_cells.sound_velocity.current_value,
                                            DataContainer().numeric.a_pseudo,
                                            DataContainer().numeric.b_pseudo)
        mock_compute_pseudo.assert_any_call(dt,
                                            self.mock_disc.additional_dof_density.current_value,
                                            self.mock_disc.additional_dof_density.new_value,
                                            self.mock_disc.right_part_size.new_value,
                                            self.mock_disc.additional_dof_sound_velocity.current_value,
                                            DataContainer().numeric.a_pseudo,
                                            DataContainer().numeric.b_pseudo)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionCell, "compute_time_step", spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_enriched_elements_new_time_step(self, mock_compute_dt, mock_disc_list):
        """
        Test de la m�thode compute_enriched_elements_new_time_step
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        mock_compute_dt.return_value = np.array([0.])

        self.my_cells.compute_enriched_elements_new_time_step()

        mock_compute_dt.assert_any_call(DataContainer().numeric.cfl,
                                        DataContainer().numeric.cfl_pseudo,
                                        self.my_cells.density.current_value,
                                        self.my_cells.density.new_value,
                                        self.mock_disc.left_part_size.new_value,
                                        self.my_cells.sound_velocity.new_value,
                                        self.my_cells.pseudo.current_value,
                                        self.my_cells.pseudo.new_value)
        mock_compute_dt.assert_called_with(DataContainer().numeric.cfl,
                                           DataContainer().numeric.cfl_pseudo,
                                           self.mock_disc.additional_dof_density.current_value,
                                           self.mock_disc.additional_dof_density.new_value,
                                           self.mock_disc.right_part_size.new_value,
                                           self.mock_disc.additional_dof_sound_velocity.new_value,
                                           self.mock_disc.additional_dof_artificial_viscosity.current_value,
                                           self.mock_disc.additional_dof_artificial_viscosity.new_value)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_stress_tensor(self, mock_disc_list):
        """
        Test de la m�thode compute_enriched_stress_tensor
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        exact_stress_gauche = np.array([[1., -1., 2.]])
        exact_stress_droite = np.array([[1., -1., 2.]])

        # Test pour la partie gauche
        self.my_cells._deviatoric_stress_new = np.array([[2., 0., 3.]])
        self.my_cells.pressure.new_value = np.array([1.])
        self.my_cells.pseudo.new_value = np.array([0.])

        # Test pour la partie droite
        self.mock_disc.additional_dof_deviatoric_stress_new = np.array([[2., 0., 3.]])
        self.mock_disc.additional_dof_pressure.new_value = np.array([1.])
        self.mock_disc.additional_dof_artificial_viscosity.new_value = np.array([0.])

        self.my_cells.compute_enriched_stress_tensor()

        # V�rification
        np.testing.assert_allclose(self.my_cells.stress, exact_stress_gauche)
        np.testing.assert_allclose(self.mock_disc.additional_dof_stress, exact_stress_droite)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionHansboEnrichedCell, "compute_discontinuity_borders_velocity",
                       spec=classmethod, new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionCell, "general_method_deviator_strain_rate",
                       spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_enriched_deviatoric_strain_rate(self, mock_compute_D, mock_disc_borders,
                                                     mock_disc_list):
        """
        Test de la m�thode compute_enriched_deviatoric_strain_rate
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        dt = 1.
        node_coord_new = np.array([[0.,], [1.,]])
        node_velocity_new = np.array([[-1, ], [1., ]])
        u_disc_g = np.array([-0.5])
        u_disc_d = np.array([0.5])
        mock_disc_borders.return_value = u_disc_g, u_disc_d
        mock_compute_D.return_value = np.array([[1., 1., 1.], ])
        self.mock_disc.plastic_cells = True

        self.my_cells.compute_enriched_deviatoric_strain_rate(dt, node_coord_new, node_velocity_new)

        mock_disc_borders.assert_called_with(self.mock_disc, node_velocity_new)
        mock_compute_D.assert_called()


    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionCell, "general_method_deviator_strain_rate", spec=classmethod,
                       new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionHansboEnrichedCell, "compute_enriched_shear_modulus",
                       spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_enriched_deviatoric_stress_tensor(self, mock_compute_G,
                                                       mock_compute_D, mock_disc_list):
        """
        Test de la m�thode compute_enriched_deviatoric_stress_tensor
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        mask = np.array([True])
        dt = 1.
        coord_noeud_new = np.array([[-1.], [0, ]])
        vitesse_noeud_new = np.array([[50, ], [-20, ]])
        # Mock topo :
        topologie = mock.MagicMock(Topology1D)
        type(topologie).cells_in_contact_with_node = mock.PropertyMock(
            return_value=np.array([[-1, 0], [0, 1], [1, 2], [2, -1]]))
        mock_compute_D.return_value = np.array([[1., -0.5, -0.5]])

        self.my_cells._deviatoric_stress_current = np.array([[0., 0., 0.]])
        self.my_cells.shear_modulus.new_value = np.array([3.])
        self.my_cells.shear_modulus.current_value = np.array([3.])
        self.my_cells._deviatoric_strain_rate = np.array([[1., -0.5, -0.5]])
        exact_new_S_left = np.array([[6., -3., -3.]])

        self.mock_disc.additional_dof_deviatoric_stress_current = np.array([[0., 0., 0.]])
        self.mock_disc.additional_dof_shear_modulus.new_value = np.array([1.])
        self.mock_disc.additional_dof_deviatoric_strain_rate = np.array([[1., -0.5, -0.5]])
        exact_new_S_right = np.array([[2., -1, -1]])

        self.my_cells.compute_enriched_deviatoric_stress_tensor(coord_noeud_new,
                                                                vitesse_noeud_new, dt)

        mock_compute_G.assert_called_with()
        mock_compute_D.assert_called()
        np.testing.assert_allclose(self.my_cells._deviatoric_stress_new, exact_new_S_left)
        np.testing.assert_allclose(self.mock_disc._additional_dof_deviatoric_stress_new,
                                   exact_new_S_right)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_shear_modulus(self, mock_disc_list):
        """
        Test de la m�thode compute_enriched_shear_modulus
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]

        self.my_cells.compute_enriched_shear_modulus()

        np.testing.assert_allclose(self.my_cells.shear_modulus.new_value,
                                   DataContainer().material_target.initial_values.shear_modulus_init)
        np.testing.assert_allclose(self.mock_disc.additional_dof_shear_modulus.new_value,
                                   DataContainer().material_target.initial_values.shear_modulus_init)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_equivalent_plastic_strain_rate(self, mock_disc_list):
        """
        Test de la m�thode compute_enriched_equivalent_plastic_strain_rate
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        mask = np.array([True])
        dt = 1.
        self.mock_disc.plastic_cells = True
        exact_equivalent_plastic_strain_rate_left = np.array([0.66666666666666663])
        exact_equivalent_plastic_strain_rate_right = np.array([0.33333333333333331])

        self.my_cells._deviatoric_stress_new = np.array([[2., -1., -1.]])
        self.my_cells.shear_modulus.current_value = np.array([1.])
        self.my_cells.yield_stress.current_value = np.array([1.])
        self.mock_disc.additional_dof_deviatoric_stress_new = np.array([[4, -2., -2.]])
        self.mock_disc.additional_dof_shear_modulus.current_value = np.array([4.])
        self.mock_disc.additional_dof_yield_stress.current_value = np.array([2.])

        self.my_cells.compute_enriched_equivalent_plastic_strain_rate(mask, dt)

        np.testing.assert_allclose(self.my_cells.equivalent_plastic_strain_rate,
                                   exact_equivalent_plastic_strain_rate_left)
        np.testing.assert_allclose(
            self.mock_disc._additional_dof_equivalent_plastic_strain_rate,
            exact_equivalent_plastic_strain_rate_right)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_apply_plastic_correction_on_enriched_deviatoric_stress_tensor(self, mock_disc_list):
        """
        Test dela m�thode apply_plastic_correction_on_enriched_deviatoric_stress_tensor
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]
        mask = np.array([True])
        self.mock_disc.plastic_cells = True
        # set partie gauche
        self.my_cells._deviatoric_stress_new = np.array([[10., -5., -5.], ])
        self.my_cells.yield_stress.current_value = np.array([10.])
        exact_deviator_plastic = np.array([[6.666667, -3.333333, -3.333333]])
        # set partie droite
        self.mock_disc._additional_dof_deviatoric_stress_new = np.array([[5., -2.5, -2.5], ])
        self.mock_disc.additional_dof_yield_stress.current_value = np.array([10.])
        exact_deviator_plastic_right = np.array([[8.006408, -4.003204, -4.003204]])

        self.my_cells.apply_plastic_correction_on_enriched_deviatoric_stress_tensor(mask)

        np.testing.assert_allclose(self.my_cells._deviatoric_stress_new,
                                   exact_deviator_plastic, rtol=1.e-5)
        np.testing.assert_allclose(self.mock_disc._additional_dof_deviatoric_stress_new,
                                   exact_deviator_plastic_right, rtol=1.e-5)

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_yield_stress(self, mock_disc_list):
        """
        Test de la m�thode compute_enriched_shear_modulus
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_disc]

        self.my_cells.compute_enriched_yield_stress()

        np.testing.assert_allclose(self.my_cells.yield_stress.new_value,
                                   DataContainer().material_target.initial_values.yield_stress_init)
        np.testing.assert_allclose(self.mock_disc.additional_dof_yield_stress.new_value,
                                   DataContainer().material_target.initial_values.yield_stress_init)
