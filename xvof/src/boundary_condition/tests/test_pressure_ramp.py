#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Unittest for PressureRamp module
"""
import unittest
from xvof.src.boundary_condition.pressure_ramp import PressureRamp


class RampPressureTest(unittest.TestCase):
    """
    Test case utilis� pour test les fonctions du module 'PressureRamp'
    """
    def setUp(self):
        """
        Tests initialization
        """
        self._first_value = 2.
        self._second_value = 4.
        self._start_time = 1.
        self._end_time = 3.
        self._pressure = PressureRamp(self._first_value, self._second_value, self._start_time, self._end_time)


    def test_evaluate_first_value(self):
        """
        Tests the evaluate method of the PressureRamp class for t < start_time
        """
        print __name__ + " : Test evaluate first value"
        value_test = self._pressure.evaluate(self._start_time - 10.)
        self.assertEqual(value_test, self._first_value)
        print "__[OK]"
        

    def test_evaluate_second_value(self):
        """
        Tests the evaluate method of the PressureRamp class for t > end_time
        """
        print __name__ + " : Test evaluate second value"
        value_test = self._pressure.evaluate(self._end_time + 10.)
        self.assertEqual(value_test, self._second_value)
        print "__[OK]"

    def test_evaluate_medium_value(self):
        """
        Tests the evaluate method of the PressureRamp class for t > (start_time + end_time) / 2 
        """
        print __name__ + " : Test evaluate medium value"
        value_test = self._pressure.evaluate(0.5 * (self._start_time + self._end_time))
        self.assertEqual(value_test, 0.5 * (self._first_value + self._second_value))
        print "__[OK]"


if __name__ == '__main__':
    unittest.main()