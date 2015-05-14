#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base abstraite (interface) définissant un critère de rupture
"""
from abc import ABCMeta, abstractmethod


class RuptureCriterion():
    """
    Une interface pour les critères de rupture
    """
    # Nécessaire pour spécifier l'interface
    __metaclass__ = ABCMeta
    #

    def __init__(self):
        pass

    @abstractmethod
    def checkCriterion(self, cell, *args, **kwargs):
        """
        Vérification du critère de rupture sur la maille
        passée en argument
        """
        pass
