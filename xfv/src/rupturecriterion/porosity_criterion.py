# -*- coding: utf-8 -*-
"""
Implementation of PorosityCriterion class
"""
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion


class PorosityCriterion(RuptureCriterion):   # pylint: disable=too-few-public-methods
    """
    A rupture criterion based on porosity value
    """
    def __init__(self, p_limit):
        super().__init__()
        self.__limit_porosity = p_limit

    def check_criterion(self, cells, *args, **kwargs):
        """
        Return the mask of the cells where porosity is above the threshold porosity

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where porosity is above the threshold porosity
        """
        return cells.porosity.new_value > self.__limit_porosity
