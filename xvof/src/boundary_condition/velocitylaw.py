#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant une loi de vitesse
"""
from abc import abstractmethod
from xvof.boundary_condition.boundary_condition import BoundaryCondition


class VelocityLaw(BoundaryCondition):
    """
    Une interface pour les lois de vitesse
    """

    def __init__(self):
        super(VelocityLaw, self).__init__()
        self._type = "velocity"

    @abstractmethod
    def evaluate(self, time, *args, **kwargs):
        """
        Evaluation de la vitesse au temps <time>
        """
        pass
