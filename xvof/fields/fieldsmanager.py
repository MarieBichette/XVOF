#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant le gestionnaire de champ
"""
from xvof.fields.enrichedfield import EnrichedField
from xvof.fields.field import Field


class FieldManager(object):
    '''
    Gestionnaire de champ
    '''
    def __init__(self):
        self.__fields = {}

    def addClassicalField(self, name, size, current_value=0., new_value=0.):
        '''
        Ajoute un champ au gestionnaire
        '''
        if name not in self.__fields.keys():
            self.__fields[name] = Field(size, current_value, new_value)
        else:
            raise KeyError('Le champ {:s} existe déjà dans le gestionnaire!'.format(name))

    def moveClassicalToEnrichedFields(self, size):
        '''
        Transforme un champ classique en un champ enrichi
        '''
        for name, field in self.__fields.items():
            self.__fields[name] = EnrichedField(size, field.current_value, field.new_value)

    def getField(self, name):
        '''
        Retourne le champ demandé
        '''
        return self.__fields[name]

    def incrementFields(self):
        '''
        Incrémente tous les champs
        '''
        for field in self.__fields.values():
            field.incrementValues()

    def printInfos(self):
        '''
        Quelques infos
        '''
        for name, field in self.__fields.items():
            print "<-- Champ {:s} de type {}-->".format(name, type(field))
