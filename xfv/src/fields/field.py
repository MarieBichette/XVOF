# -*- coding: iso-8859-1 -*-
"""
Implementing the Field class
"""
import numpy as np


class Field:
    """
    Classical physical field on cells or nodes. Owning current and future values in numpy.arrays
    """

    def __init__(self, size, current_value=0., new_value=0.):
        """
        :param size: size of arrays (i.e nodes or cells number)
        :param current_value: current field value
        :param new_value: future field value

        :type size: int
        :type current_value: float or numpy.array
        :type new_value: float or numpy.array
        """
        self.__size = size
        self.__current = np.empty([size], dtype=np.float64, order='C')
        self.__future = np.empty([size], dtype=np.float64, order='C')
        self.__current[:] = current_value
        self.__future[:] = new_value

    def __str__(self):
        """
        :return: informations about the field
        """
        return "{:s} of size {:d}".format(self.__class__.__name__, self.size)

    @property
    def size(self):
        """
        :return: the size of the field (i.e number of cells or nodes on which the field is defined)
        """
        return self.__size

    @property
    def current_value(self):
        """
        :return: a copy of current field value
        :rtype: numpy.array
        """
        return self.__current[:]

    @current_value.setter
    def current_value(self, value):
        """
        Set value as the current value of the field (useful for Hansbo enrichment)

        :param value: new field value to set
        :type value: float or numpy.array
        """
        self.__current[:] = value

    @property
    def new_value(self):
        """
        :return: a copy of the future field value
        :rtype: numpy.array
        """
        return self.__future[:]

    @new_value.setter
    def new_value(self, value):
        """
        Set value as the future value of the field

        :param value: new field value to set
        :type value: float or numpy.array
        """
        self.__future[:] = value

    def increment_values(self):
        """
        Increment field values
        """
        self.__current[:] = self.__future[:]
