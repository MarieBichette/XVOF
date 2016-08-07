# -*- coding: iso-8859-1 -*-
"""
Implementing the EnrichedField class
"""
import numpy as np

from xvof.fields.field import Field


def from_geometry_to_enrich_field(left_field, right_field):
    """
    :return: Enriched field
    :rtype: numpy.array
    """
    return (right_field - left_field) * 0.5


def from_geometry_to_classic_field(left_field, right_field):
    """
    :return: Classical field
    :rtype: numpy.array
    """
    return (right_field + left_field) * 0.5


def from_enrich_to_left_part_field(classic_field, enriched_field):
    """
    :return: Left field
    :rtype: numpy.array
    """
    return classic_field - enriched_field


def from_enrich_to_right_part_field(classic_field, enriched_field):
    """
    :return: Right field
    :rtype: numpy.array
    """
    return classic_field + enriched_field


class EnrichedField(Field):
    """
    Physical field = classical physical field + enriched field
    """
    def __init__(self, size, current_value, new_value):
        super(EnrichedField, self).__init__(size, current_value, new_value)
        self.__enr_current = np.empty([size], dtype=np.float64, order='C')
        self.__enr_future = np.empty([size], dtype=np.float64, order='C')
        self.__enr_current[:] = 0.
        self.__enr_future[:] = 0.

    def incrementValues(self):
        """
        Increment field values
        """
        super(EnrichedField, self).incrementValues()
        self.__enr_current[:] = self.__enr_future[:]

    @property
    def current_enr_value(self):
        """
        :return: a copy of the enriched current field values
        :rtype: numpy.array
        """
        return self.__enr_current[:]

    @property
    def new_enr_value(self):
        """
        :return: a copy of the enriched future field values
        :rtype: numpy.array
        """
        return self.__enr_future[:]

    @new_enr_value.setter
    def new_enr_value(self, value):
        """
        Set value as the future value of the enriched field

        :param value: new value to set
        :type value: float ou numpy.array
        """
        self.__enr_future[:] = value

    @property
    def current_left_value(self):
        """
        :return: Current left field
        """
        return from_enrich_to_left_part_field(self.current_value, self.current_enr_value)

    @property
    def current_right_value(self):
        """
        :return: Current right field
        """
        return from_enrich_to_right_part_field(self.current_value, self.current_enr_value)

    @property
    def new_left_value(self):
        """
        :return: Future left field
        """
        return from_enrich_to_left_part_field(self.new_value, self.new_enr_value)

    @property
    def new_right_value(self):
        """
        :return: Future right field
        """
        return from_enrich_to_right_part_field(self.new_value, self.new_enr_value)
