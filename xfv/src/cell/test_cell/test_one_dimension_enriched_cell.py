# -*- coding: utf-8 -*-
"""
Classe de test du module element1dupgraded
"""
import unittest
import unittest.mock as mock
import numpy as np
import os
from xfv.src.cell.one_dimension_enriched_cell import OneDimensionEnrichedCell
from xfv.src.data.data_container import DataContainer
from xfv.src.discontinuity.discontinuity import Discontinuity


class OneDimensionEnrichedCellTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Tests setup for class
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_epp.json")
        DataContainer(data_file_path)

    @classmethod
    def tearDownClass(cls):
        DataContainer.clear()
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

    def setUp(self):
        """
        Test preparation
        """
        self.nbr_cells = 4
        self.my_elements = OneDimensionEnrichedCell(self.nbr_cells)

    def tearDown(self):
        pass

    def test_classical(self):
        """
        Test of the property lassical
        """
        np.testing.assert_array_equal(self.my_elements.classical,
                                      np.ones(self.nbr_cells, dtype="bool"))
        self.my_elements.classical[0] = False
        np.testing.assert_array_equal(self.my_elements.classical,
                                      np.array([False, True, True, True]))
        self.my_elements.classical[0] = True

    def test_enriched(self):
        """
        Test of the property enriched
        """
        np.testing.assert_array_equal(self.my_elements.enriched,
                                      np.zeros(self.nbr_cells, dtype="bool"))
        self.my_elements.classical[0] = False
        np.testing.assert_array_equal(self.my_elements.enriched,
                                      np.array([True, False, False, False]))
        self.my_elements.classical[0] = True

    def test_compute_new_left_right_size(self):
        """
        Test of compute_new_left_right_size
        """
        disc = mock.MagicMock(Discontinuity)
        type(disc.left_part_size).current_value = mock.PropertyMock(return_value=np.array([1.0]))
        type(disc.right_part_size).current_value = mock.PropertyMock(return_value=np.array([1.0]))
        u1h = 0.
        u2h = 2.
        ug = -0.5
        ud = 1.0
        OneDimensionEnrichedCell.compute_new_left_right_size(0.5, disc, u1h, u2h, ug, ud)
        np.testing.assert_array_almost_equal(disc.left_part_size.new_value, np.array([0.75]))
        np.testing.assert_array_almost_equal(disc.right_part_size.new_value, np.array([1.5]))

    def test_compute_new_left_right_density(self):
        """
        Test of compute_new_left_right_density
        """
        disc = mock.MagicMock(Discontinuity)
        type(disc.left_part_size).current_value = mock.PropertyMock(return_value=np.array([1.]))
        type(disc.left_part_size).new_value = mock.PropertyMock(return_value=np.array([0.5]))
        type(disc.right_part_size).current_value = mock.PropertyMock(return_value=np.array([1.]))
        type(disc.right_part_size).new_value = mock.PropertyMock(return_value=np.array([2.]))
        density_left = 1.
        density_right = 1.
        density_left_new, density_right_new = OneDimensionEnrichedCell.\
            compute_new_left_right_density(density_left, density_right, disc)
        np.testing.assert_array_almost_equal(density_left_new, np.array([2.]))
        np.testing.assert_array_almost_equal(density_right_new, np.array([0.5]))

if __name__ == "__main__":
    unittest.main()
