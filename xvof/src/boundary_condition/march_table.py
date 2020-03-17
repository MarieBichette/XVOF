#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une pression à partir d'un fichier texte en entrée (table de marche)
"""
import numpy as np
from xvof.src.boundary_condition.pressurelaw import PressureLaw


class MarchTablePressure(PressureLaw):
    """
    Une pression à deux paliers constants
    """
    def __init__(self, data_file):
        super(MarchTablePressure, self).__init__()
        self.__data_file = data_file
        self.__data = np.loadtxt(self.__data_file)
        print "Pressure is imposed from file" + self.__data_file
        self.time_data = self.__data[:, 0]
        self.pressure_data = self.__data[:, 1]

    def evaluate(self, time, *args, **kwargs):
        """
        Evalue la pression à retourner en fonction du temps par rapport à la table de marche
        :param time: temps courant
        :return: la pression à imposer
        """
        if time > self.time_data[-1]:
            pressure = 0
        else:
            indice_limite_inf = np.where(self.time_data <= time)[-1][0]
            indice_limite_sup = np.where(self.time_data >= time)[0][0]
            # import ipdb ; ipdb.set_trace()
            if indice_limite_inf == indice_limite_sup:
                # le temps courant se situe dans la table de marche
                pressure = self.pressure_data[indice_limite_inf]
            else:
                # on interpole linéairement la pression entre les deux temps présents dans la table de marche
                pente = (self.pressure_data[indice_limite_sup] - self.pressure_data[indice_limite_inf]) / \
                        (self.time_data[indice_limite_sup] - self.time_data[indice_limite_inf])
                pressure = self.pressure_data[indice_limite_inf] + (time - self.time_data[indice_limite_inf]) * pente
        return pressure

