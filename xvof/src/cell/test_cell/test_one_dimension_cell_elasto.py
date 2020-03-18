#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import numpy as np
import unittest
import mock
import os
from xvof.src.cell.cell import Cell
from xvof.src.cell.one_dimension_cell import OneDimensionCell
from xvof.src.mesh.topology1d import Topology1D
from xvof.src.data.data_container import DataContainer
from xvof.src.rheology.constantshearmodulus import ConstantShearModulus
from xvof.src.cell.test_cell.test_variables import TestVariables

class OneDimensionCellTest(unittest.TestCase):

    def setUp(self):
        data_file_path = os.path.realpath(os.path.join(os.getcwd(), "../tests/0_UNITTEST/XDATA_elasto.xml"))
        self.test_datacontainer = DataContainer(data_file_path)
        self.nbr_cells = 4
        self.my_cells = OneDimensionCell(self.nbr_cells)
        self.my_cells.cell_in_target = np.ones(self.nbr_cells, dtype='bool')

        self.test_variables = TestVariables(4, 5)
        self.test_variables.define_elasto_variables()

        self.mask = np.array([True, True, False, False])

    def tearDown(self):
        pass

    @mock.patch.object(OneDimensionCell, "apply_equation_of_state", spec=classmethod, new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionCell, "add_elastic_energy_method", spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_new_pressure_with_elasticity(self, mock_add_elasticity, mock_apply_eos):
        """
        Test de la méthode compute_new_pressure
        """
        print "Test " + __name__
        # L'option élasticité est activée dans le jeu de donnée
        type(DataContainer().material_target.constitutive_model).elasticity_model = \
            mock.PropertyMock(return_value=ConstantShearModulus)

        # Configuration des mocks
        mask_classic = self.mask
        mock_apply_eos.return_value[mask_classic] = self.my_cells.energy.new_value[mask_classic], \
                                      self.my_cells.pressure.new_value[mask_classic], \
                                      self.my_cells.sound_velocity.new_value[mask_classic]
        mock_add_elasticity.return_value[mask_classic] = self.my_cells.energy.new_value[mask_classic]
        # les valeurs de return permettent juste à la méthode de ne pas planter mais ne sont pas importantes...
        dt = 1.

        self.my_cells.compute_new_pressure(mask_classic, dt)

        mock_add_elasticity.assert_called_with(dt,
                                               self.my_cells.density.current_value[mask_classic],
                                               self.my_cells.density.new_value[mask_classic],
                                               mock.ANY,
                                               mock.ANY,
                                               mock.ANY)
        # /!\ mock.ANY sert à éviter de faire appel aux quantités déviatoriques car la méthode assert_called_with ne
        # fonctionne pas sur les arguments de type array de taille > 1 (ValueError). On vérifie quand même que la
        # méthode add_elasticity est appelée avec une partie des arguments correcte

        mock_apply_eos.assert_called_with(self.my_cells,
                                          DataContainer().material_target.constitutive_model.eos,
                                          self.my_cells.density.current_value[mask_classic],
                                          self.my_cells.density.new_value[mask_classic],
                                          self.my_cells.pressure.current_value[mask_classic],
                                          self.my_cells.pressure.new_value[mask_classic],
                                          self.my_cells.energy.current_value[mask_classic],
                                          self.my_cells.energy.new_value[mask_classic],
                                          self.my_cells.pseudo.current_value[mask_classic],
                                          self.my_cells.sound_velocity.new_value[mask_classic])
        print "C'est la meilleure vérification qu'on puisse faire mais pas sure que ça suffise..."
        print "__[OK]"

    def test_compute_shear_modulus(self):
        """
        Test de la méthode compute_shear_modulus
        """
        print "Test " + __name__
        self.my_cells.compute_shear_modulus()
        np.testing.assert_allclose(self.my_cells.shear_modulus.new_value, self.test_variables.shear_modulus_init)
        print "__[OK]"

    def test_compute_yield_stress(self):
        """
        Test de la méthode compute_yield_stress
        """
        print "Test " + __name__
        self.my_cells.compute_yield_stress()
        np.testing.assert_allclose(self.my_cells.yield_stress.new_value, self.test_variables.yield_stress_init)
        print "__[OK]"

    def test_compute_complete_stress_tensor(self):
        """
        Test de la méthode compute_complete_stress_tensor
        """
        print "Test " + __name__
        mask = self.mask
        self.my_cells._stress = np.copy(self.test_variables.stress_old)
        self.my_cells._deviatoric_stress_new = np.copy(self.test_variables.deviatoric_stress_new)
        self.my_cells.pressure.new_value = np.copy(self.test_variables.pressure_new)
        self.my_cells.pseudo.new_value = np.copy(self.test_variables.pseudo_new)
        self.my_cells.compute_complete_stress_tensor(mask)
        np.testing.assert_allclose(self.my_cells.stress[mask], self.test_variables.stress_new[mask])
        np.testing.assert_allclose(self.my_cells.stress[~mask], self.test_variables.stress_old[~mask])
        print "__[OK]"

    @mock.patch.object(OneDimensionCell, "compute_shear_modulus", spec=classmethod, new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionCell, "compute_deviator_strain_rate", spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_deviatoric_stress_tensor(self, mock_compute_D, mock_compute_G):
        """
        Test de la méthode compute_deviatoric_stress_tensor
        """
        print __name__ + " : Test compute deviatoric stress"
        mask = self.mask
        dt = self.test_variables.dt
        coord_noeud_new = np.copy(self.test_variables.node_coord_new)
        vitesse_noeud_new = np.copy(self.test_variables.node_velocity_new)
        # Mock topo :
        topologie = mock.MagicMock(Topology1D)
        type(topologie).cells_in_contact_with_node = mock.PropertyMock(
            return_value=np.array([[-1, 0], [0, 1], [1, 2], [2, 3], [3, -1]]))

        self.my_cells._deviatoric_stress_current = np.copy(self.test_variables.deviatoric_stress_old)
        self.my_cells.shear_modulus.current_value = np.copy(self.test_variables.shear_modulus_old)
        self.my_cells._deviatoric_strain_rate = np.copy(self.test_variables.strain_rate_dev_new)
        mock_compute_D.return_value = self.my_cells._deviatoric_strain_rate

        self.my_cells.compute_deviatoric_stress_tensor(mask, topologie, coord_noeud_new, vitesse_noeud_new, dt)

        mock_compute_G.assert_called_with()
        mock_compute_D.assert_called_with(mask, dt, topologie, coord_noeud_new, vitesse_noeud_new)
        np.testing.assert_allclose(self.my_cells._deviatoric_stress_new[mask],
                                   self.test_variables.deviatoric_stress_new[mask])
        np.testing.assert_allclose(self.my_cells._deviatoric_stress_new[~mask],
                                   self.test_variables.deviatoric_stress_old[~mask])
        print "__[OK]"

    @mock.patch.object(Topology1D, "nodes_belonging_to_cell", new_callable=mock.PropertyMock,
                       return_value=np.array([[0, 1], [1, 2], [2, 3], [3, 4]]))
    def test_compute_deviator_strain_rate(self, mock_topologie):
        """
        Test de la méthode  compute_deviatoric_strain_rate
        """
        print __name__ + " : Test compute deviator strain rate"
        # données d'entrée
        mask = self.mask
        dt = self.test_variables.dt
        coord_noeud_new = np.copy(self.test_variables.node_coord_new)
        vitesse_noeud_new = np.copy(self.test_variables.node_velocity_new)
        topo_ex = Topology1D(5, 4)
        # Test de la méthode compute_deviator
        self.my_cells.compute_deviator_strain_rate(mask, dt, topo_ex, coord_noeud_new, vitesse_noeud_new)
        np.testing.assert_allclose(self.my_cells._deviatoric_strain_rate[mask],
                                   self.test_variables.strain_rate_dev_new[mask])
        print "__[OK]"

if __name__ == "__main__":
    unittest.main()
