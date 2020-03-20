#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import numpy as np
import unittest
import mock
import os
from xvof.src.cell.one_dimension_cell import OneDimensionCell
from xvof.src.mesh.topology1d import Topology1D
from xvof.src.data.data_container import DataContainer


class OneDimensionCellTest(unittest.TestCase):

    def setUp(self):
        data_file_path = os.path.join(os.path.dirname(__file__), "../../../tests/0_UNITTEST/XDATA_hydro.xml")
        self.test_datacontainer = DataContainer(data_file_path)
        self.nbr_cells = 4
        self.my_cells = OneDimensionCell(self.nbr_cells)
        self.my_cells.cell_in_target = np.ones(self.nbr_cells, dtype='bool')

    def tearDown(self):
        DataContainer.clear()
        pass

    def test_compute_size(self):
        """
        Test of compute_size method
        """
        topo = mock.MagicMock(Topology1D)
        type(topo).nodes_belonging_to_cell = mock.PropertyMock(return_value=np.array([[0, 1], [1, 2], [2, 3], [3, 4]]))
        node_coord = np.array([[-0.35, ], [0.15, ], [0.2, ], [0.25, ], [0.5, ]])
        self.my_cells.compute_size(topo, node_coord)
        np.testing.assert_allclose(self.my_cells.size_t, np.array([0.5, 0.05, 0.05, 0.25]))

    def test_compute_new_size(self):
        """
        Test of compute_new_size method
        """
        topo = mock.MagicMock(Topology1D)
        type(topo).nodes_belonging_to_cell = mock.PropertyMock(return_value=np.array([[0, 1], [1, 2], [2, 3],[3, 4]]))
        mask = np.array([True, True, False, False])
        self.my_cells._size_t = np.array([0.5, 0.05, 0.05, 0.25])
        node_new_coord = np.array([[-0.25, ], [0.1, ], [0.2, ], [0.45, ], [0.85, ]])
        self.my_cells.compute_new_size(topo, node_new_coord, mask)
        np.testing.assert_allclose(self.my_cells.size_t_plus_dt, np.array([0.35, 0.1, 0., 0.]))

    def test_compute_mass(self):
        """
        Test of compute_mass method
        """
        assert self.test_datacontainer.geometric.section == 0.000003141592653589793
        self.my_cells.density.current_value = np.array([8129., 8129., 8129., 8129.])
        self.my_cells._size_t = np.array([0.6, 0.1, 0.15, 0.6])
        self.my_cells.compute_mass()
        np.testing.assert_allclose(self.my_cells.mass,
                                   np.array([0.015323,  0.002554,  0.003831,  0.015323]), rtol=1.e-4)

    def test_compute_new_density(self):
        """
        Test of compute_new_density method
        """
        self.my_cells.density.current_value = np.ones(self.nbr_cells) * 8930.
        self.my_cells._size_t = np.array([0.6, 0.1, 0.15, 0.6])
        self.my_cells._size_t_plus_dt = np.array([0.8, 0.1, 0.27, 0.8])
        mask = np.array([True, True, True, False])
        self.my_cells.compute_new_density(mask)
        np.testing.assert_allclose(self.my_cells.density.new_value, np.array([6697.5,  8930., 4961.111111, 8930.]))

    def test_compute_new_pressure_without_elasticity(self):
        """
        Test de la méthode compute_new_pressure
        """
        mask = np.array([True, True, True, False])
        dt = 1.

        self.my_cells.density.current_value = np.ones(self.nbr_cells) * 8930.
        self.my_cells.pressure.current_value = np.ones(self.nbr_cells) * 1.e+5
        self.my_cells.energy.current_value = np.ones(self.nbr_cells) * 6.5

        self.my_cells.density.new_value = np.array([9000., 8900., 8915., 8920.])
        self.my_cells.sound_velocity.new_value = np.zeros([self.nbr_cells])
        self.my_cells.pressure.new_value = np.zeros([self.nbr_cells])
        self.my_cells.energy.new_value = np.zeros([self.nbr_cells])
        self.my_cells.pseudo.new_value = np.zeros([self.nbr_cells])

        self.my_cells.compute_new_pressure(mask, dt)
        np.testing.assert_allclose(self.my_cells.energy.new_value, np.array([487.203942, 94.056026, 28.379115, 0.]))
        np.testing.assert_allclose(self.my_cells.pressure.new_value, np.array([1.103734e+09, -4.640127e+08,
                                                                               -2.323423e+08, 0.]), rtol=1e-5)

    def test_compute_new_pseudo(self):
        """
        Test of compute_new_pseudo
        """
        mask = np.array([False, True, True, True])
        self.my_cells.density.current_value = np.array([8700., 3200, 2171, 8700.])
        self.my_cells.density.new_value = np.array([8800., 3500, 2175, 8500.])
        self.my_cells.sound_velocity.current_value = np.array([4400, 3200, 1140, 4400])
        self.my_cells._size_t_plus_dt = np.array([0.025, 0.01, 0.005, 0.025])
        self.my_cells.pseudo.new_value = np.ones(self.nbr_cells) * -1
        assert self.test_datacontainer.numeric.a_pseudo == 1.2
        assert self.test_datacontainer.numeric.b_pseudo == 0.25
        self.my_cells.compute_new_pseudo(1.2e-08, mask)
        np.testing.assert_allclose(self.my_cells.pseudo.new_value,
                                   np.array([-1., 2.25427729e+13, 2.00897590e+09, 0.]))

    def test_compute_new_time_step(self):
        """
        Test de la méthode compute_enriched_elements_new_time_step
        """
        self.my_cells.density.current_value = np.array([8000., 8500., 9500., 8000.])
        self.my_cells.density.new_value = np.array([8120., 8440., 9000., 8120.])
        self.my_cells.sound_velocity.new_value = np.array([100., 100., 100., 100.])
        self.my_cells.pseudo.current_value = np.array([0., 0., 0., 0.])
        self.my_cells.pseudo.new_value = np.array([1., 1., 1., 1.])
        self.my_cells._size_t_plus_dt = np.array([1., 1., 2., 1.])
        assert self.test_datacontainer.numeric.cfl == 0.34
        assert self.test_datacontainer.numeric.cfl_pseudo == 0.
        self.my_cells._dt = np.array([1., 2., 5., 6.])
        mask = np.array([False, True, True, False])
        self.my_cells.compute_new_time_step(mask)
        np.testing.assert_allclose(self.my_cells._dt, np.array([1., 3.400000e-03, 6.800000e-03, 6.]))

    def test_compute_complete_stress_tensor_hydro(self):
        """
        Test de la méthode compute_complete_stress_tensor
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
        expected_result = np.array([[-999., -999, -999.],
                                    [-23, -23, -23],
                                    [1000., -500., -500.],
                                    [4000., -2000., -2000.]])
        np.testing.assert_allclose(self.my_cells.stress, expected_result)

    def test_impose_pressure(self):
        """
        Test de impose_pressure
        """
        self.my_cells.pressure.new_value = np.array([2000., 4000, 8000., 4000.])
        self.my_cells.impose_pressure(0, -1.)
        np.testing.assert_array_equal(self.my_cells.pressure.new_value, np.array([-1, 4000, 8000., 4000.]))

if __name__ == "__main__":
    unittest.main()
