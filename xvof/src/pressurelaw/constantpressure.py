#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une pression constante
"""
from xvof.src.pressurelaw.pressurelaw import PressureLaw


class ConstantPressure(PressureLaw):
    """
    Une pression constante
    """
    def __init__(self, value):
        super(ConstantPressure, self).__init__()
        self.__value = value

    def evaluate(self, time, *args, **kwargs):
        return self.__value
