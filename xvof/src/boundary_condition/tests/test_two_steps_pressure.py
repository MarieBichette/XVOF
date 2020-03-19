# -*- coding: utf-8 -*-
"""
Unittests of the two_steps_pressure module
"""
import unittest

import xvof.src.boundary_condition.two_steps_pressure as TwoStpPressure

class TwoStepsPressureTest(unittest.TestCase):
    """
    This class tests the methods of the TwoStepsPressure class
    """
    def setUp(self):
        """
        Tests initialization
        """
        self._first_value = 2.
        self._second_value = 4.
        self._critical_time = 0.
        self._pressure = TwoStpPressure.TwoStepsPressure(self._first_value,
                                                         self._second_value,
                                                         self._critical_time)

    def test_evaluate_fisrt_value(self):
        """
        Tests the evaluate method for t < critical_time
        """
        value_test = self._pressure.evaluate(self._critical_time - 10.)
        self.assertEqual(value_test, self._first_value)


    def test_evaluate_second_value(self):
        """
        Tests the evaluate method for t > critical_time
        """
        value_test = self._pressure.evaluate(self._critical_time + 10.)
        self.assertEqual(value_test, self._second_value)


if __name__ == '__main__':
    unittest.main()
