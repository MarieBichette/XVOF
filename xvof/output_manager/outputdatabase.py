"""
Implementing the Hdf5Database class giving access to hdf5 storage
"""

import h5py
import numpy as np


class OutputDatabase(object):
    """
    A class to store simulation fields in an hdf5 database
    """

    def __init__(self, path_to_hdf5):
        self.__db = h5py.File(path_to_hdf5, 'w')
        self.__current_group = None
        self.__nb_sav = 0

    def add_time(self, time):
        """
        Create an hdf5 group containing all fields for the current time
        """
        self.__db.flush()  # Flushing to print preceding time steps
        self.__current_group = self.__db.create_group("{:9.8g}".format(time))
        self.__nb_sav += 1

    def add_field(self, field_name, values, **kwargs):
        """
        Create a dataset corresponding to field_name and storing the values. All extra keywords arguments
        are stored as attributes of the dataset
        """
        ds = self.__current_group.create_dataset(field_name, values.shape, np.result_type(values))
        ds[...] = values
        for k, v in kwargs.items():
            ds.attrs[k] = v

    def close(self):
        """
        Close the database
        """
        self.__db.flush()
        self.__db.close()
