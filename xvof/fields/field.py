# -*- coding: iso-8859-1 -*-
"""
Implementing the Field class
"""
import numpy as np


class Field(object):
    """
    Classical physical field on cells or nodes. Owning current and future values in numpy.arrays
    """
    def __init__(self, size, current_value, new_value):
        """
        :param size: size of arrays (i.e nodes or cells number)
        :param current_value: current field value
        :param new_value: future field value

        :type size: int
        :type current_value: float or numpy.array
        :type new_value: float or numpy.array
        """
        self.__values = {'current': np.empty([size], dtype=np.float64, order='C'),
                         'new': np.empty([size], dtype=np.float64, order='C')}
        self.__values['current'][:] = current_value
        self.__values['new'][:] = new_value

    def incrementValues(self):
        """
        Increment field values
        """
        self.__values['current'][:] = self.__values['new'][:]

    @property
    def current_value(self):
        """
        :return: a copy of current field value
        :rtype: numpy.array
        """
        return self.__values['current'][:]

    @property
    def new_value(self):
        """
        :return: a copy of the future field value
        :rtype: numpy.array
        """
        return self.__values['new'][:]

    @new_value.setter
    def new_value(self, value):
        """
        Set value as the future value of the field

        :param value: new field value to set
        :type value: float or numpy.array
        """
        self.__values['new'][:] = value
