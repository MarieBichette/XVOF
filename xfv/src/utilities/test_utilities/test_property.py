#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Classe de test des property dans les differentes classes
Valable pour toutes les classes qui ont des property
"""
import unittest
import numpy as np


class PropertyTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions property de renvoi
    """
    def setUp(self):
        """
        Préparation des tests unitaires
        """
        class MyClass:
            """
            A fake class with some properties
            """
            def __init__(self):
                """
                Construction d'une classe quelconque contenant 2 proprétés _prop et _prop_enriched
                """
                self._field = np.array([0., 1., 2.])
                self._field_enriched = np.array([0., 4., 0.])

            def set_field(self, value):
                """
                Sets the self._field to the imposed value
                :param value: value to impose
                """
                self._field = value

            def set_enriched_field(self, value):
                """
                Sets the self._field_enriched to the imposed value
                :param value: value to impose
                """
                self._field_enriched = value

            @property
            def prop(self):
                """
                une méthode qui renvoie _prop (simple renvoi)
                """
                return self._field

            @property
            def prop_copy(self):
                """
                une méthode qui renvoie _prop (simple renvoi)
                """
                return np.copy(self._field)

            @property
            def prop_total(self):
                """
                une méthode qui renvoie _prop (simple renvoi) et une méthode qui renvoie le champ
                total (calcul dans prop) - 1ere ecriture
                """
                res = self._field
                res += self._field_enriched
                return res

            @property
            def prop_total_point(self):
                """
                une méthode qui renvoie _prop (simple renvoi) et une méthode qui renvoie le champ
                total (calcul dans prop) - 2eme ecriture
                """
                res = self._field[:]
                res += self._field_enriched
                return res

            @property
            def prop_total_copy(self):
                """
                une méthode qui renvoie _prop (simple renvoi) et une méthode qui renvoiele champ
                total (calcul dans prop) - 3eme ecriture
                """
                res = np.copy(self._field)
                res += self._field_enriched
                return res

        self.my_class = MyClass()

    def test_prop(self):
        """
        Teste la méthode de renvoi simple
        """
        # Remise à zéro de la classe
        self.my_class.set_field(np.array([0., 1., 2.]))
        self.my_class.set_enriched_field(np.array([0., 4., 0.]))
        # test de ...
        np.testing.assert_array_equal(self.my_class.prop,
                                      self.my_class._field)  # pylint: disable=protected-access
        result = self.my_class.prop
        result[0] += 10
        np.testing.assert_array_equal(result, np.array([10., 1., 2.]))
        np.testing.assert_array_equal(self.my_class._field,  # pylint: disable=protected-access
                                      np.array([10., 1., 2.]))

    def test_prop_copy(self):
        """
        Teste la méthode de renvoi simple - ecriture copy
        """
        # Remise à zéro de la classe
        self.my_class.set_field(np.array([0., 1., 2.]))
        self.my_class.set_enriched_field(np.array([0., 4., 0.]))
        # test de ...
        np.testing.assert_array_equal(self.my_class.prop_copy,
                                      self.my_class._field)  # pylint: disable=protected-access
        result = self.my_class.prop_copy
        result[0] += 10
        np.testing.assert_array_equal(result, np.array([10., 1., 2.]))
        np.testing.assert_array_equal(self.my_class._field,  # pylint: disable=protected-access
                                      np.array([0., 1., 2.]))

    def test_prop_total(self):
        """
        Teste la méthode de renvoi avec operation - ecriture simple
        """
        # Remise à zéro de la classe
        self.my_class.set_field(np.array([0., 1., 2.]))
        self.my_class.set_enriched_field(np.array([0., 4., 0.]))
        # test de ...
        result = self.my_class.prop_total
        np.testing.assert_array_equal(result, np.array([0., 5., 2.]))
        np.testing.assert_array_equal(self.my_class._field,  # pylint: disable=protected-access
                                      np.array([0., 5., 2.]))

    def test_prop_total_point(self):
        """
        Teste la méthode de renvoi avec operation - ecriture [:]
        """
        # Remise à zéro de la classe
        self.my_class.set_field(np.array([0., 1., 2.]))
        self.my_class.set_enriched_field(np.array([0., 4., 0.]))
        # test de ...
        result = self.my_class.prop_total_point
        np.testing.assert_array_equal(result, np.array([0., 5., 2.]))
        np.testing.assert_array_equal(self.my_class._field,  # pylint: disable=protected-access
                                      np.array([0., 5., 2.]))

    def test_prop_total_copy(self):
        """
        Teste la méthode de renvoi avec operation - ecriture np.copy()
        """
        # Remise à zéro de la classe
        self.my_class.set_field(np.array([0., 1., 2.]))
        self.my_class.set_enriched_field(np.array([0., 4., 0.]))
        # test de ...
        result = self.my_class.prop_total_copy
        np.testing.assert_array_equal(result, np.array([0., 5., 2.]))
        np.testing.assert_array_equal(self.my_class._field,  # pylint: disable=protected-access
                                      np.array([0., 1., 2.]))


if __name__ == '__main__':
    unittest.main()
