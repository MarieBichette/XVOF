# -*- coding: utf-8 -*-
"""
Unittest for pressure_ramp module
"""
import unittest
from xvof.src.custom_functions.ramp import Ramp


class RampTest(unittest.TestCase):
    """
    Test the Ramp class
    """
    def setUp(self):
        """
        Tests initialization
        """
        self._first_value = 2.
        self._second_value = 4.
        self._start_time = 1.
        self._end_time = 3.
        self._pressure = Ramp(self._first_value, self._second_value,
                              self._start_time, self._end_time)

    def test_evaluate_first_value(self):
        """
        Tests the evaluate method of the Ramp class for t < start_time
        """
        value_test = self._pressure.evaluate(self._start_time - 10.)
        self.assertEqual(value_test, self._first_value)

    def test_evaluate_second_value(self):
        """
        Tests the evaluate method of the Ramp class for t > end_time
        """
        value_test = self._pressure.evaluate(self._end_time + 10.)
        self.assertEqual(value_test, self._second_value)

    def test_evaluate_medium_value(self):
        """
        Tests the evaluate method of the Ramp class for t > (start_time + end_time) / 2
        """
        value_test = self._pressure.evaluate(0.5 * (self._start_time + self._end_time))
        self.assertEqual(value_test, 0.5 * (self._first_value + self._second_value))


if __name__ == '__main__':
    unittest.main()
