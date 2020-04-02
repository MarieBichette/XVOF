# -*- coding: utf-8 -*-
# pylint: disable=protected-access
"""
one_dimension_cell module unit tests
"""
import unittest
import numpy as np
import os
from xfv.src.cell.one_dimension_cell import OneDimensionCell
from xfv.src.mesh.topology1d import Topology1D
from xfv.src.data.data_container import DataContainer


class OneDimensionCellElastoTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Tests setup for class
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_elasto.json")
        DataContainer(data_file_path)

    @classmethod
    def tearDownClass(cls):
        DataContainer.clear()
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

    def setUp(self):
        self.nbr_cells = 4
        self.my_cells = OneDimensionCell(self.nbr_cells)
        self.my_cells.cell_in_target = np.ones(self.nbr_cells, dtype='bool')
        self.mask = np.array([True, True, False, False])

    def tearDown(self):
        pass

    def test_compute_new_pressure_with_elasticity(self):
        """
        Test de la m�thode compute_new_pressure
        """
        mask_classic = np.array([True, True, False, False])
        self.my_cells.density.current_value = np.ones(self.nbr_cells) * 8930.
        self.my_cells.pressure.current_value = np.ones(self.nbr_cells) * 1.e+5
        self.my_cells.energy.current_value = np.ones(self.nbr_cells) * 7.5

        self.my_cells._deviatoric_stress_new = np.array([[200., -100, -100.],
                                                         [5., -2.5, -2.5],
                                                         [100., -50., -50.],
                                                         [400., -200., -200.]])

        self.my_cells._deviatoric_strain_rate = \
            np.array([[-0.235294,  0.117647,  0.117647],
                      [0.4444444, -0.2222222, -0.2222222],
                      [0.,  0.,  0.],
                      [0.,  0.,  0.]])

        self.my_cells.density.new_value = np.array([8940., 8970, 8920., 9000.])
        self.my_cells.sound_velocity.new_value = np.zeros([self.nbr_cells])
        self.my_cells.pressure.new_value = np.zeros([self.nbr_cells])
        self.my_cells.energy.new_value = np.zeros([self.nbr_cells])
        self.my_cells.pseudo.new_value = np.zeros([self.nbr_cells])
        delta_t = 1.

        self.my_cells.compute_new_pressure(mask_classic, delta_t)
        np.testing.assert_allclose(self.my_cells.pressure.new_value,
                                   np.array([1.557157e+08, 6.266033e+08,
                                             0.000000e+00, 0.000000e+00]), rtol=1.e-5)
        np.testing.assert_allclose(self.my_cells.energy.new_value,
                                   np.array([17.254756,  163.9763, 0., 0.]))
        np.testing.assert_allclose(self.my_cells.sound_velocity.new_value,
                                   np.array([3948.726929, 3974.84139, 0., 0.]))

    def test_compute_shear_modulus(self):
        """
        Test de la m�thode compute_shear_modulus
        """
        self.my_cells.compute_shear_modulus()
        expected_value = DataContainer().material_target.initial_values.shear_modulus_init
        np.testing.assert_allclose(self.my_cells.shear_modulus.new_value,
                                   np.ones([self.nbr_cells]) * expected_value)

    def test_compute_yield_stress(self):
        """
        Test de la m�thode compute_yield_stress
        """
        self.my_cells.compute_yield_stress()
        expected_value = DataContainer().material_target.initial_values.yield_stress_init
        np.testing.assert_allclose(self.my_cells.yield_stress.new_value,
                                   np.ones([self.nbr_cells]) * expected_value)

    def test_compute_complete_stress_tensor(self):
        """
        Test de la m�thode compute_complete_stress_tensor
        """
        mask = np.array([True, True, False, False])
        self.my_cells._stress = np.array([[2000., -1000, -1000.],
                                          [50., -25., -25.],
                                          [1000., -500., -500.],
                                          [4000., -2000., -2000.]])
        self.my_cells._deviatoric_stress_new = np.array([[200., -100, -100.],
                                                         [5., -2.5, -2.5],
                                                         [100., -50., -50.],
                                                         [400., -200., -200.]])
        self.my_cells.pressure.new_value = np.array([1000, 25, 500, 2000])
        self.my_cells.pseudo.new_value = np.array([-1, -2, -3, -4])
        self.my_cells.compute_complete_stress_tensor(mask)
        expected_result = np.array([[-799., -1099, -1099.],
                                     [-18, -25.5, -25.5],
                                     [1000., -500., -500.],
                                     [4000., -2000., -2000.]])
        np.testing.assert_allclose(self.my_cells.stress, expected_result)

    def test_compute_deviatoric_stress_tensor(self):
        """
        Test de la m�thode compute_deviatoric_stress_tensor
        """
        mask = np.array([True, True, True, False])
        delta_t = 1.
        coord_noeud_new = np.array([[-0.25, ], [0.1, ], [0.2, ], [0.45, ], [0.85, ]])
        vitesse_noeud_new = np.array([[0.1, ], [-0.05, ], [0., ], [0.2, ], [0.3, ]])
        topo_ex = Topology1D(5, 4)
        self.my_cells._deviatoric_stress_current = np.array([[2000., -1000, -1000.],
                                                             [50., -25., -25.],
                                                             [1000., -500., -500.],
                                                             [4000., -2000., -2000.]])
        self.my_cells.shear_modulus.current_value = np.array([2., 4., 6., 8.])
        self.my_cells.compute_deviatoric_stress_tensor(mask, topo_ex, coord_noeud_new,
                                                       vitesse_noeud_new, delta_t)
        np.testing.assert_allclose(self.my_cells._deviatoric_stress_new,
                                   np.array([[1999.058824,  -999.529412,  -999.529412],
                                             [53.555556,   -26.777778,   -26.777778],
                                             [1010.66666667,  -505.33333333,  -505.33333333],
                                             [0., 0., 0.]]), rtol=1.e-05)

    def test_compute_deviator_strain_rate(self):
        """
        Test de la m�thode  compute_deviatoric_strain_rate
        """
        # donn�es d'entr�e
        mask = np.array([True, True, False, False])
        delta_t = 1.
        coord_noeud_new = np.array([[-0.25, ], [0.1, ], [0.2, ], [0.45, ], [0.85, ]])
        vitesse_noeud_new = np.array([[0.1, ], [-0.05, ], [0., ], [0.2, ], [0.3, ]])
        topo_ex = Topology1D(5, 4)
        # Test de la m�thode compute_deviator
        self.my_cells.compute_deviator_strain_rate(mask, delta_t, topo_ex, coord_noeud_new,
                                                   vitesse_noeud_new)
        np.testing.assert_allclose(self.my_cells._deviatoric_strain_rate,
                                   np.array([[-0.235294,  0.117647,  0.117647],
                                             [0.4444444, -0.2222222, -0.2222222],
                                             [0., 0., 0.], [0., 0., 0.]]), rtol=1.e-05)


if __name__ == "__main__":
    unittest.main()
