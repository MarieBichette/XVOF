#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Unit tests for enriched_field class
"""
import unittest
import numpy as np
from xfv.src.fields.field import Field


class FieldTest(unittest.TestCase):
    """
    Classe test pour Field.
    """

    def setUp(self):
        """
        Tests setup
        """
        self.my_field = Field(1, np.array([0.]), np.array([0.]))

    def test_increment_values(self):
        """
        Test of the increment_values method
        """
        self.my_field.new_value = np.array([5.])
        self.my_field.increment_values()
        np.testing.assert_array_equal(self.my_field.current_value, self.my_field.new_value)

    def test_new_enr_value_setter(self):
        """
        Test of the new_value attribute setter
        """
        self.my_field.new_value = np.array([10.])
        np.testing.assert_array_equal(self.my_field.new_value, np.array([10.]))


if __name__ == "__main__":
    unittest.main()

