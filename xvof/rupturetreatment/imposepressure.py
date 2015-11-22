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
        super(ImposePressure, self).__init__()
        self.__imposed_pressure = pressure

    def applyTreatment(self, cells, ruptured_cells, *args, **kwargs):
        cells.imposePressure(ruptured_cells, self.__imposed_pressure)
