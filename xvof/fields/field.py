#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un champ
"""
import numpy as np


class Field(object):
    '''
    Champ physique classique. Il est définit par une valeur actuelle et une valeur
    future. Il s'agit en fait de deux tableaux numpy.
    '''
    def __init__(self, size, current_value, new_value):
        """
        :param size: taille des tableaux
        :param current_value: valeur actuelle à fixer
        :param new_value: valeur future à fixer

        :type size: int
        :type current_value: float ou numpy.array
        :type new_value: float ou numpy.array
        """
        self.__values = {'current': np.empty([size], dtype=np.float64, order='C'),
                         'new': np.empty([size], dtype=np.float64, order='C')}
        if isinstance(current_value, float):
            self.__values['current'].fill(current_value)
        else:
            self.__values['current'][:] = current_value[:] 
        if isinstance(new_value, float):
            self.__values['new'].fill(new_value)
        else:
            self.__values['new'][:] = new_value[:] 

    def incrementValues(self):
        '''
        Incrémente les valeurs du champ
        '''
        self.__values['current'][:] = self.__values['new'][:]

    @property
    def current_value(self):
        '''
        Valeur actuelle (au temps t) du champ

        :return: une copie des valeurs du champs au temps courant
        :rtype: numpy.array
        '''
        return self.__values['current'][:]

    @property
    def new_value(self):
        '''
        Nouvelle valeur (au temps t+dt) du champ

        :return: une copie des valeurs du champs au temps futur
        :rtype: numpy.array
        '''
        return self.__values['new'][:]

    @new_value.setter
    def new_value(self, value):
        '''
        Setter de la nouvelle valeur du champ

        :param value: nouvelle valeur du champ à fixer
        :type value: float ou numpy.array
        '''
        if isinstance(value, float):
            self.__values['new'].fill(value)
        else:
            self.__values['new'][:] = value[:]
