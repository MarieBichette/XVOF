#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module equationofstate
"""
import unittest

import xvof.equationsofstate.equationofstate as eos


class EquationOfStateTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'EquationOfState'
    """
    def test_init(self):
        """ Test du constructeur EquationOfState() """
        with self.assertRaises(TypeError):
            eos.EquationOfState()

if __name__ == '__main__':
    unittest.main()
