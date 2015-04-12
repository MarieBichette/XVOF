#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module equationofstate
"""
import xvof.equationsofstate.equationofstate as eos
import unittest


class EquationOfStateTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'EquationOfState'
    """
    def test_init(self):
        """
        Test du constructeur
        """
        with self.assertRaises(TypeError):
            eos.EquationOfState()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
if __name__ == '__main__':
    unittest.main()