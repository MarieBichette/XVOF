#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test des property dans les differentes classes
Valable pour toutes les classes qui ont des property
"""
import numpy as np
import unittest

class PropertyTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions property de renvoi
    """
    def setUp(self):
        """
        Préparation des tests unitaires
        """
        class MyClass(object):
            def __init__(self):
                """construction d'une classe quelconque contenant 2 proprétés _prop et _prop_enriched """
                self._field = np.array([0., 1., 2.])
                self._field_enriched = np.array([0., 4., 0.])

            @property
            def prop(self):
                """une méthode qui renvoie _prop (simple renvoi) """
                return self._field

            @property
            def prop_copy(self):
                """une méthode qui renvoie _prop (simple renvoi) """
                return np.copy(self._field)

            @property
            def prop_total(self):
                """une méthode qui renvoie _prop (simple renvoi) et une méthode qui renvoiele champ total
                (calcul dans prop) - 1ere ecriture"""
                res = self._field
                res += self._field_enriched
                return res

            @property
            def prop_total_point(self):
                """une méthode qui renvoie _prop (simple renvoi) et une méthode qui renvoiele champ total
                (calcul dans prop) - 2eme ecriture"""
                res = self._field[:]
                res += self._field_enriched
                return res

            @property
            def prop_total_copy(self):
                """une méthode qui renvoie _prop (simple renvoi) et une méthode qui renvoiele champ total
                (calcul dans prop) - 3eme ecriture"""
                res = np.copy(self._field)
                res += self._field_enriched
                return res

        self.my_class = MyClass()

    def test_prop(self):
        """Teste la méthode de renvoi simple"""
        np.testing.assert_array_equal(self.my_class.prop, self.my_class._field)
        self.my_class.prop[0] += 10.
        np.testing.assert_array_equal(self.my_class.prop, np.array([10., 1., 2.]))
        np.testing.assert_array_equal(self.my_class._field, np.array([0., 1., 2.]))

    def test_prop_copy(self):
        """Teste la méthode de renvoi simple - ecriture copy"""
        np.testing.assert_array_equal(self.my_class.prop_copy, self.my_class._field)
        self.my_class.prop[0] += 10.
        np.testing.assert_array_equal(self.my_class.prop, np.array([10., 1., 2.]))
        np.testing.assert_array_equal(self.my_class._field, np.array([0., 1., 2.]))

    def test_prop_total(self):
        """Teste la méthode de renvoi avec operation - ecriture simple"""
        np.testing.assert_array_equal(self.my_class.prop_total, np.array([0., 5., 2.]))
        appel1 = self.my_class.prop_total # appel dela prop quelque part dans le code
        np.testing.assert_array_equal(self.my_class._field, np.array([0., 1., 2.]))

    def test_prop_total_point(self):
        """Teste la méthode de renvoi avec operation - ecriture [:]"""
        np.testing.assert_array_equal(self.my_class.prop_total_point, np.array([0., 5., 2.]))
        appel1 = self.my_class.prop_total_point # appel dela prop quelque part dans le code
        np.testing.assert_array_equal(self.my_class._field, np.array([0., 1., 2.]))

    def test_prop_total_copy(self):
        """Teste la méthode de renvoi avec operation - ecriture np.copy()"""
        np.testing.assert_array_equal(self.my_class.prop_total_copy, np.array([0., 5., 2.]))
        appel1 = self.my_class.prop_total_copy # appel dela prop quelque part dans le code
        np.testing.assert_array_equal(self.my_class._field, np.array([0., 1., 2.]))



if __name__ == '__main__':
    unittest.main()