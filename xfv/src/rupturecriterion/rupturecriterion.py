#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementing the RuptureCriterion abstract base class
"""
from abc import ABCMeta, abstractmethod


class RuptureCriterion(object, metaclass=ABCMeta):  # pylint: disable=too-few-public-methods
    """
    An abstract base class for rupture criteria
    """

    def __init__(self):
        pass

    @abstractmethod
    def check_criterion(self, cells, *args, **kwargs):
        """
        Check of the rupture criterion on the cells in arguments
        :param cells: cells on which to check the criterion
        """
        pass
