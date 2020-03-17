# -*- coding: iso-8859-1 -*-
"""
Implementing the EnrichedField class
"""
import numpy as np

from xvof.src.fields.field import Field


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
        self.__enr_current = np.zeros([size], dtype=np.float64, order='C')
        self.__enr_future = np.zeros([size], dtype=np.float64, order='C')

    def increment_values(self):
        """
        Increment field values
        """
        super(EnrichedField, self).increment_values()
        self.__enr_current[:] = self.__enr_future[:]
