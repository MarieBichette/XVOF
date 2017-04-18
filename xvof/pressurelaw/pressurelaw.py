#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base abstraite (interface) définissant une loi de pression
"""
from abc import ABCMeta, abstractmethod


class PressureLaw(object):
    """
    Une interface pour les lois de pression
    """
    # Nécessaire pour spécifier l'interface
    __metaclass__ = ABCMeta
    #

    def __init__(self):
        pass

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):
        """
        Evaluation de la pression au temps <time>
        """
        pass
