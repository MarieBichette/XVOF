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


# def write_evolution_from_db(path_to_hdf5_db, output_file_name, item, id_item, field):
#     """
#     Save the hdf5 database in a .dat file
#     :param path_to_hdf5_db : path where the hdf5 database is stored
#     :param item: cell or node
#     :param id_item: int, id of the item to be saved
#     :param field: item fields to be saved
#     :param output_file_name : file name for outputs
#     :return:
#     """
#     directory = os.path.dirname(os.path.abspath(path_to_hdf5_db))
#     path_to_save_data = os.path.join(directory, output_file_name)
#     print("writing file {:} with fields {} ...".format(os.path.relpath(path_to_save_data), field.title))
#     to_write = get_field_evolution_in_time_for_item(path_to_hdf5_db, id_item, field.title)
#     header = "-----Time History for {:} {:} :".format(item, str(id_item)) + os.linesep
#     header += "{:^15s}    {:^15s} ".format("Time[s]", field.label)
#     my_format = ['%+10.9e', '%+10.9e']
#
#     with open(path_to_save_data, 'w') as f:
#         np.savetxt(f, to_write, fmt=my_format, header=header)
#         f.close()
#
#
# def write_profile_from_db(path_to_hdf5_db, output_file_name, time, field):
#     """
#     Save the hdf5 database in a .dat file
#     :param path_to_hdf5_db : path where the hdf5 database is stored
#     :param output_file_name : file to write
#     :param time : time
#     :param field: fields to be saved
#     :return:
#     """
#     directory = os.path.dirname(os.path.abspath(path_to_hdf5_db))
#     path_to_save_data = os.path.join(directory, output_file_name)
#     my_hd = OutputDatabaseExploit(path_to_hdf5_db)
#     array_to_write = my_hd.extract_true_field_at_time(field.title, time)[:]
#
#     np.savetxt(path_to_save_data, array_to_write, fmt=['%+9.8e','%+9.8e'])
#     print("writing file {:} with fields {} ...".format(path_to_save_data, field.title))





