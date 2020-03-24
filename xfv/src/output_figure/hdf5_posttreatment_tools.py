#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-
"""
Tools for posttreat the hdf5 band to write the true fields for cells and nodes selected
"""

import os
import numpy as np
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit


def read_database_at_time_t(path_to_hdf5_db, id_item, field, t):
    my_hd = OutputDatabaseExploit(path_to_hdf5_db)
    fields_item_history_at_time = np.zeros([2])
    fields_item_history_at_time[0] = t
    if "Classical" in field:
        fields_item_history_at_time[1] = my_hd.extract_field_at_time(field, t)[id_item]
    elif "Additional" in field:
        fields_item_history_at_time[1] = my_hd.extract_field_at_time(field, t)[id_item][1]
    else:
        fields_item_history_at_time[1] = my_hd.extract_true_field_at_time(field, t)[id_item][1]
    return fields_item_history_at_time


def read_database_in_array(path_to_hdf5_db, id_item, field, modulo=1):
    """
    Read the hdf5 database in path_to_hdf5 and returns an array with time and field in field_list
    :param path_to_hdf5_db: str, path to the hdf5 database
    :param id_item: int, id of the items to be posttreated
    :param field: list of Time_field, fields to be posttreated. (à rentrer comme dans time_figure_tools...)
    :return: fields_item_history : np.array(time, fields in field list)
    """
    my_hd = OutputDatabaseExploit(path_to_hdf5_db)
    fields_item_history = []
    index_time = 0

    # Au début, on remplit normalement, pour tous les temps
    for i in range(0, len(my_hd.saved_times) / 4):
    # for t in my_hd.saved_times:
        t = my_hd.saved_times[i]
        fields_item_history.append(read_database_at_time_t(path_to_hdf5_db, id_item, field, t).tolist())
        # fields_item_history[index_time, :] =
        index_time += 1

    # A la fin, on remplit une valeur sur modulo
    for i in range(len(my_hd.saved_times) / 4, len(my_hd.saved_times)):
        if i % modulo == 0:
            t = my_hd.saved_times[i]
            # fields_item_history[index_time, :] = read_database_at_time_t(path_to_hdf5_db, id_item, field, t)
            fields_item_history.append(read_database_at_time_t(path_to_hdf5_db, id_item, field, t).tolist())
            index_time += 1
    return np.array(fields_item_history)


def write_evolution_from_db(path_to_hdf5_db, output_file_name, item, id_item, field):
    """
    Save the hdf5 database in a .dat file
    :param path_to_hdf5_db : path where the hdf5 database is stored
    :param item: cell or node
    :param id_item: int, id of the item to be saved
    :param field: item fields to be saved
    :param output_file_name : file name for outputs
    :return:
    """
    directory = os.path.dirname(os.path.abspath(path_to_hdf5_db))
    path_to_save_data = os.path.join(directory, output_file_name)
    print("writing file {:} with fields {} ...".format(os.path.relpath(path_to_save_data), field.title))
    to_write = read_database_in_array(path_to_hdf5_db, id_item, field.title)
    header = "-----Time History for {:} {:} :".format(item, str(id_item)) + os.linesep
    header += "{:^15s}    {:^15s} ".format("Time[s]", field.label)
    my_format = ['%+10.9e', '%+10.9e']

    with open(path_to_save_data, 'w') as f:
        np.savetxt(f, to_write, fmt=my_format, header=header)
        f.close()


def write_profile_from_db(path_to_hdf5_db, output_file_name, time, field):
    """
    Save the hdf5 database in a .dat file
    :param path_to_hdf5_db : path where the hdf5 database is stored
    :param output_file_name : file to write
    :param time : time
    :param field: fields to be saved
    :return:
    """
    directory = os.path.dirname(os.path.abspath(path_to_hdf5_db))
    path_to_save_data = os.path.join(directory, output_file_name)
    my_hd = OutputDatabaseExploit(path_to_hdf5_db)
    array_to_write = my_hd.extract_true_field_at_time(field.title, time)[:]

    np.savetxt(path_to_save_data, array_to_write, fmt=['%+9.8e','%+9.8e'])
    print("writing file {:} with fields {} ...".format(path_to_save_data, field.title))





