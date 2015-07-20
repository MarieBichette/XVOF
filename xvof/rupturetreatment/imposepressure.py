#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un traitement de la rupture par imposition de la pression
"""
from xvof.rupturetreatment.rupturetreatment import RuptureTreatment


class ImposePressure(RuptureTreatment):
    """
    Un traitement de rupture qui impose la pression
    """
    def __init__(self, pressure):
        RuptureTreatment.__init__(self)
        self.__imposed_pressure = pressure

    def applyTreatment(self, cell, *args, **kwargs):
        cell.imposePressure(self.__imposed_pressure)
