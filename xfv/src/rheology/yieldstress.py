# -*- coding: utf-8 -*-
"""
Implementation d'une classe de limite d'élasticité (interface)
"""


class YieldStress:  # pylint: disable=too-few-public-methods
    """
    Interface for yield stress computation
    """
    def __init__(self, initial_value):
        """
        Initialization of the class
        :param initial_value: yield_stress initial value
        """
        self.yield_stress = initial_value

    @classmethod
    def compute(cls) -> float:
        """
        Compute the new value of shear modulus
        :return : float
        """
