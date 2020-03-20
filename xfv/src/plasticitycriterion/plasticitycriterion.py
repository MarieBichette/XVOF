#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementing the PlasticityCriterion abstract base class
"""
from abc import ABCMeta, abstractmethod


class PlasticityCriterion(object):
    """
    An abstract base class for plasticity criteria
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def checkCriterion(self, cells, *args, **kwargs):
        """
        Check of the plasticity criterion on the cells in arguments

        :param cells: cells on which to check the criterion
        """
        pass
