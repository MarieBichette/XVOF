#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une condition limite générique
"""
from abc import abstractmethod


class BoundaryCondition(object):
    """
    Une condition limite générique
    """
    def __init__(self):
        self._type = None

    def type(self):
        """
        Accesseur sur le type de condition limite (pressure ou velocity)
        :return:
        """
        return self._type

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):
        pass
