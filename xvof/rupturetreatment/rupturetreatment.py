#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base abstraite (interface) définissant un traitement de rupture
"""
from abc import ABCMeta, abstractmethod


class RuptureTreatment(object):
    """
    Une interface pour les traitements de la rupture
    """
    # Nécessaire pour spécifier l'interface
    __metaclass__ = ABCMeta
    #

    def __init__(self):
        pass

    @abstractmethod
    def applyTreatment(self, cell, *args, **kwargs):
        """
        Application du traitement de rupture sur la maille
        passée en argument
        """
        pass
