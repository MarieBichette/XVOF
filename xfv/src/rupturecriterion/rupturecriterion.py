#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementing the RuptureCriterion abstract base class
"""
from abc import ABCMeta, abstractmethod


class RuptureCriterion(object):  # pylint: disable=too-few-public-methods
    """
    An abstract base class for rupture criteria
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def check_criterion(self, cells, *args, **kwargs):
        """
        Check of the rupture criterion on the cells in arguments
        :param cells: cells on which to check the criterion
        """
        pass
