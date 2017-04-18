#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Unit tests for enriched_field class
"""
import numpy as np
import unittest
import mock


from xvof.fields.field import Field


class FieldTest(unittest.TestCase):
    """
    Classe test pour Field.
    Les proprétés sont testées dans utilities
    """

    def setUp(self):
        """
        Tests setup
        """
        self.my_field = Field(1, np.array([0.]), np.array([0.]))


    def test_incrementValues(self):
        self.my_field.new_value = np.array([5.]) # ou mettre en mock ???
        self.my_field.incrementValues()
        np.testing.assert_array_equal(self.my_field.current_value, self.my_field.new_value)

    def test_new_enr_value_setter(self):
        self.my_field.new_value(np.array([10.]))
        np.testing.assert_array_equal(self.my_field.new_value, np.array([10.]))



if __name__ == "__main__":
    unittest.main()

