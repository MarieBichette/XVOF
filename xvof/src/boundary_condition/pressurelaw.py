#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant une loi de pression
"""
from abc import abstractmethod
from xvof.src.boundary_condition.boundary_condition import BoundaryCondition


class PressureLaw(BoundaryCondition):
    """
    Une interface pour les lois de pression
    """

    def __init__(self):
        super(PressureLaw, self).__init__()
        self._type = "pressure"

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):
        """
        Evaluation de la pression au temps <time>
        """
        pass
