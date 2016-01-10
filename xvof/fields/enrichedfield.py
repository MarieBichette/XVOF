#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un champ enrichi comme
composition de deux champs classiques
"""
from xvof.fields.field import Field
import numpy as np


class EnrichedField(Field):
    '''
    Champ physique = champ physique classique + enrichi
    '''
    def __init__(self, size, current_value, new_value):
        super(EnrichedField, self).__init__(size, current_value, new_value)
        self.__enr_values = {'current': np.empty([size], dtype=np.float64, order='C'),
                         'new': np.empty([size], dtype=np.float64, order='C')}
        self.__enr_values['current'][:] = 0.
        self.__enr_values['new'][:] = 0.

    def incrementValues(self):
        '''
        Incrémente les valeurs du champ
        '''
        super(EnrichedField, self).incrementValues()
        self.__enr_values['current'][:] = self.__enr_values['new'][:]
        
        
    @property
    def current_enr_value(self):
        '''
        Valeur actuelle (au temps t) du champ

        :return: une copie des valeurs du champs au temps courant
        :rtype: numpy.array
        '''
        return self.__enr_values['current'][:]

    @property
    def new_enr_value(self):
        '''
        Nouvelle valeur (au temps t+dt) du champ

        :return: une copie des valeurs du champs au temps futur
        :rtype: numpy.array
        '''
        return self.__enr_values['new'][:]

    @new_enr_value.setter
    def new_enr_value(self, value):
        '''
        Setter de la nouvelle valeur du champ

        :param value: nouvelle valeur du champ à fixer
        :type value: float ou numpy.array
        '''
        self.__enr_values['new'][:] = value

    @classmethod
    def fromGeometryToEnrichField(cls, champ_gauche, champ_droite):
        """
        Renvoi le champ enrichi à partir des champs gauche et droite
        """
        return (champ_droite - champ_gauche) * 0.5

    @classmethod
    def fromGeometryToClassicField(cls, champ_gauche, champ_droite):
        """
        Renvoi le champ classique à partir des champs gauche et droite
        """
        return (champ_droite + champ_gauche) * 0.5

    @classmethod
    def fromEnrichToLeftPartField(cls, champ_classic, champ_enrich):
        """
        Renvoi le champ à gauche d'après les champs classsique et enrichis
        """
        return champ_classic - champ_enrich

    @classmethod
    def fromEnrichToRightPartField(cls, champ_classic, champ_enrich):
        """
        Renvoi le champ à droite d'après les champs classsique et enrichis
        """
        return champ_classic + champ_enrich

    @property
    def current_left_value(self):
        '''
        Renvoie la valeur courante du champ à gauche dans l'élément enrichi
        '''
        return EnrichedField.fromEnrichToLeftPartField(self.current_value, self.current_enr_value)

    @property
    def current_right_value(self):
        '''
        Renvoie la valeur courante du champ à droite dans l'élément enrichi
        '''
        return EnrichedField.fromEnrichToRightPartField(self.current_value, self.current_enr_value)

    @property
    def new_left_value(self):
        '''
        Renvoie la nouvelle valeur du champ à gauche dans l'élément enrichi
        '''
        return EnrichedField.fromEnrichToLeftPartField(self.new_value, self.new_enr_value)

    @property
    def new_right_value(self):
        '''
        Renvoie la nouvelle valeur du  champ à droite dans l'élément enrichi
        '''
        return EnrichedField.fromEnrichToRightPartField(self.new_value, self.new_enr_value)


