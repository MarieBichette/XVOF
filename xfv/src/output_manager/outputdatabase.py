"""
Implementing the Hdf5Database class giving access to hdf5 storage
"""

import h5py


class OutputDatabase:
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

        :param time: time to be added
        """
        self.__db.flush()  # Flushing to print preceding time steps
        self.__current_group = self.__db.create_group("{:g}".format(time))
        self.__nb_sav += 1

    def add_field(self, f_name, values, **kwargs):
        """
        Create a dataset corresponding to field_name and storing the values.
        All extra keywords arguments are stored as attributes of the dataset

        :param f_name: name of the field to be added
        """
        data_set = self.__current_group.create_dataset(f_name, data=values)
        for key, value in list(kwargs.items()):
            data_set.attrs[key] = value

    def close(self):
        """
        Close the database
        """
        self.__db.flush()
        self.__db.close()
