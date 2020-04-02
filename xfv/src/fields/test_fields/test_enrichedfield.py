#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Unit tests for enriched_field class
"""
import unittest
import numpy as np
import xfv.src.fields.enrichedfield as enr_field


class EnrichedFieldTest(unittest.TestCase):
    """
    Classe test pour EnrichedField.
    """

    def setUp(self):
        """
        Tests setup
        """
        self.left_field = np.array([4.])
        self.right_field = np.array([2.])
        self.classic_field = np.array([3.])
        self.enriched_field = np.array([-1.])
        self.my_enrichedField = enr_field.EnrichedField(1, np.array([0.]), np.array([0.]))

    def test_from_geometry_to_enrich_field(self):
        """
        Test of the from_geometry_to_enrich_field translation
        """
        my_res = enr_field.from_geometry_to_enrich_field(self.left_field, self.right_field)
        np.testing.assert_array_equal(my_res, self.enriched_field)

    def test_from_geometry_to_classic_field(self):
        """
        Test of the from_geometry_to_classic_field translation
        """
        my_res = enr_field.from_geometry_to_classic_field(self.left_field, self.right_field)
        np.testing.assert_array_equal(my_res, self.classic_field)

    def test_from_enrich_to_left_part_field(self):
        """
        Test of the from_enrich_to_left_part_field translation
        """
        my_res = enr_field.from_enrich_to_left_part_field(self.classic_field, self.enriched_field)
        np.testing.assert_array_equal(my_res, self.left_field)

    def test_from_enrich_to_right_part_field(self):
        """
        Test of the from_enrich_to_right_part_field translation
        :return:
        """
        my_res = enr_field.from_enrich_to_right_part_field(self.classic_field, self.enriched_field)
        np.testing.assert_array_equal(my_res, self.right_field)


if __name__ == "__main__":
    unittest.main()
