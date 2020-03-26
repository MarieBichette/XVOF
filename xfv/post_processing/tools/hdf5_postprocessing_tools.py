#!/usr/bin/env python3.7
#  -*- coding: utf-8 -*-
"""
Tools for post processing of the hdf5 band to write the true fields for cells and nodes selected
"""

import os
import numpy as np
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit


def get_field_evolution_in_time_for_item(path_to_hdf5_db, id_item, field, modulo=1):
    """
    Read the hdf5 database in path_to_hdf5 and returns an array with time and field in field_list

    :param path_to_hdf5_db: path to the hdf5 database
    :param id_item: id of the items to be post processed
    :param field: field to be post processed
    :param modulo: get one time over modulo to increase this method performance

    :type str
    :type int
    :type string
    :type int

    :return: np.array(time, fields in field list)
    """
    my_hd = OutputDatabaseExploit(path_to_hdf5_db)
    field_item_history = np.zeros([len(my_hd.saved_times), 2])
    index_time = 0
    for i in range(len(my_hd.saved_times)):
        # if i correspond to one over modulo time of the database, get field at time t
        if i % modulo == 0:
            t = my_hd.saved_times[i]
            field_item_history[index_time, 0] = t
            field_item_history[index_time, 1] = _field_at_time_at_item(my_hd, id_item, field, t)
            index_time += 1
    return field_item_history


def get_field_profile_at_time(hdf5_band, field, t):
    """
    Extracts a profile of field at time t

    :param hdf5_band: hdf5 band containing the output results
    :param field: field to extract
    :param t: time
    :type OutputDatabaseExploit
    :type str
    :type float

    :return: tuple(coordinates, field profile)
    """
    data = hdf5_band.extract_true_field_at_time(field, t)
    coord = data[:, 0]
    field_value = data[:, 1]
    return coord, field_value


def _field_at_time_at_item(hd_band, id_item, field, t):
    """
    Extracts the field at a given item and at a given time

    :param hd_band: hdf5 band containing the output results
    :param id_item: item id to look at
    :param field: field to look at
    :param t: time to look at

    :type OutputDatabaseExploit
    :type int
    :type str
    :type float

    :return: float
    """
    if "Classical" in field:
        result = hd_band.extract_field_at_time(field, t)[id_item]
    elif "Additional" in field:
        result = hd_band.extract_field_at_time(field, t)[id_item][1]
    else:
        result = hd_band.extract_true_field_at_time(field, t)[id_item][1]
    return result
