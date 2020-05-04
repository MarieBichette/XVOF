# -*- coding: utf-8 -*-
"""
Implementing the PlasticityCriterion abstract base class
"""
from abc import ABCMeta, abstractmethod


class PlasticityCriterion(metaclass=ABCMeta):
    """
    An abstract base class for plasticity criteria
    """

    @staticmethod
    @abstractmethod
    def check_criterion(cells):
        """
        Check the plasticity criterion on the cells in arguments

        :param cells: cells on which to check the criterion
        """

    @staticmethod
    @abstractmethod
    def check_criterion_on_right_part_cells(cells):
        """
        Check the plasticity criterion on the discontinuity in arguments

        :param cells: cells on which to check the criterion
        """
