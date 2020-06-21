# -*- coding: utf-8 -*-
# pylint: disable=protected-access
"""
Test of the PenaltyContact
"""
import unittest
import numpy as np
from xfv.src.contact.penalty import PenaltyContact
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.data.enriched_mass_matrix_props import LumpMenouillardMassMatrixProps


class PenaltyContactTest(unittest.TestCase):
    """
    Test case used to test the 'PenaltyContact' module
    """
    def setUp(self):
        """
        Tests initialisation
        """
        self.test_penalty_contact = PenaltyContact(10.)
        self.disc = Discontinuity(0, np.array([True, False]), np.array([False, True]), 0.5,
                                  LumpMenouillardMassMatrixProps())
        self.disc.mass_matrix_enriched.compute_enriched_mass_matrix_left_part(0., 2., 0.5)
        self.disc.mass_matrix_enriched.compute_enriched_mass_matrix_right_part(2., 0., 0.5)

    def test_compute_penalty_force(self):
        """
        Test of the method _compute_penalty_force
        """
        result_negative = self.test_penalty_contact._compute_penalty_force(2.)
        self.assertEqual(result_negative, 0.)

        result_positive = self.test_penalty_contact._compute_penalty_force(-5.)
        self.assertEqual(result_positive, -50.)

    def test_compute_contact_force(self):
        """
        Test of the method compute_contact_force
        """
        delta_t = 1.
        node_velocity = np.array([[100., ], [50., ]])
        self.disc.discontinuity_opening.new_value = - 5.0
        contact_force = self.test_penalty_contact.compute_contact_force(node_velocity, self.disc,
                                                                        delta_t)
        self.assertEqual(contact_force, -50.0)


if __name__ == '__main__':
    unittest.main()
