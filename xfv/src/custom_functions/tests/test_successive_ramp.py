# -*- coding: utf-8 -*-
"""
Unittests of the successive_pressure_ramp module
"""
import unittest

from xfv.src.custom_functions.successive_ramp import SuccessiveRamp
from xfv.src.custom_functions.ramp import Ramp


class SuccessiveRampTest(unittest.TestCase):
    """
    Tests the SuccessiveRamp class
    """
    def setUp(self):
        """
        Tests initialization
        """
        self._first_value = 2.
        self._second_value = 4.
        self._third_value = 0.
        self._start_time = 1.
        self._end_first_slope_time = 2.
        self._end_creneau = 3.
        self._end_time = 4.
        fst_ramp = Ramp(self._first_value, self._second_value, self._start_time,
                        self._end_first_slope_time)
        snd_ramp = Ramp(self._second_value, self._third_value, self._end_creneau, self._end_time)
        self._function = SuccessiveRamp(fst_ramp, snd_ramp)

    def test_evaluate_first_value(self):
        """
        Tests the evaluate method for t < start_time
        """
        value_test = self._function.evaluate(self._start_time - 10.)
        self.assertEqual(value_test, self._first_value)

    def test_evaluate_third_value(self):
        """
        Tests the evaluate method for t > end_time
        """
        value_test = self._function.evaluate(self._end_time + 10.)
        self.assertEqual(value_test, self._third_value)

    def test_evaluate_medium_value_first(self):
        """
        Tests the evaluate method for t in the first slope
        """
        value_test = self._function.evaluate(0.5 * (self._start_time + self._end_first_slope_time))
        self.assertEqual(value_test, 0.5 * (self._first_value + self._second_value))

    def test_evaluate_plateau_value(self):
        """
        Tests the evaluate method for t in the constant part between the slopes
        """
        value_test = self._function.evaluate((self._end_first_slope_time + self._end_creneau) * 0.5)
        self.assertEqual(value_test, self._second_value)

    def test_evaluate_medium_value_second(self):
        """
        Tests the evaluate method for t in the second slope
        """
        value_test = self._function.evaluate(0.5 * (self._end_creneau + self._end_time))
        self.assertEqual(value_test, 0.5 * (self._second_value + self._third_value))


if __name__ == '__main__':
    unittest.main()
