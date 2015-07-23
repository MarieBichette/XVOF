#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un champ
"""


class Field(object):
    '''
    Champ physique classique
    '''
    def __init__(self, current, new):
        self.__values = {'current': current,
                         'new': new}

    def incrementValues(self):
        '''
        Incrémente les valeurs du champ
        '''
        self.__values['current'] = self.__values['new']

    @property
    def current_value(self):
        '''
        Valeur actuelle (au temps t) du champ
        '''
        return self.__values['current']

    @property
    def new_value(self):
        '''
        Nouvelle valeur (au temps t+dt) du champ
        '''
        return self.__values['new']

    @new_value.setter
    def new_value(self, value):
        '''
        Setter de la nouvelle valeur du champ
        '''
        self.__values['new'] = value
