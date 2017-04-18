#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Unit tests for enriched_field class
"""
import numpy as np
import unittest
import mock

import xvof.fields.enrichedfield as enr_field
from xvof.fields.field import Field


class EnrichedFieldTest(unittest.TestCase):
    """
    Classe test pour EnrichedField.
    Les proprétés sont testées dans utilities
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
        my_res = enr_field.from_geometry_to_enrich_field(self.left_field, self.right_field)
        np.testing.assert_array_equal(my_res, self.enriched_field)

    def test_from_geometry_to_classic_field(self):
        my_res = enr_field.from_geometry_to_classic_field(self.left_field, self.right_field)
        np.testing.assert_array_equal(my_res, self.classic_field)

    def test_from_enrich_to_left_part_field(self):
        my_res = enr_field.from_enrich_to_left_part_field(self.classic_field, self.enriched_field)
        np.testing.assert_array_equal(my_res, self.left_field)

    def test_from_enrich_to_right_part_field(self):
        my_res = enr_field.from_enrich_to_right_part_field(self.classic_field, self.enriched_field)
        np.testing.assert_array_equal(my_res, self.right_field)

    def test_incrementValues(self):
        self.my_enrichedField.new_enr_value = np.array([5.]) # ou mettre en mock ???
        self.my_enrichedField.incrementValues()
        np.testing.assert_array_equal(self.my_enrichedField.current_enr_value, self.my_enrichedField.new_enr_value)

    def test_new_enr_value_setter(self):
        self.my_enrichedField.new_enr_value(np.array([10.]))
        np.testing.assert_array_equal(self.my_enrichedField.new_enr_value, np.array([10.]))



if __name__ == "__main__":
    unittest.main()

