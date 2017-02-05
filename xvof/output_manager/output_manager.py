"""
Implementing the OutputManager class
"""
from collections import namedtuple
import os

from xvof.output_manager.hdf5database import Hdf5Database
from xvof.utilities.singleton import Singleton


DatabaseBuildInfos = namedtuple("DatabaseBuildInfos", ["database_object", "fields", "time_controler"])
Field = namedtuple("Field", ["name", "owner", "attr_name"])

class OutputTimeControler(object):
    """
    Compute the time or iteration of outputs
    """
    def __init__(self, identifier, time_period=None, iteration_period=None):
        self.__id = identifier
        if time_period is None and iteration_period is None:
            iteration_period = 1
        elif time_period is not None and iteration_period is not None:
            raise ValueError("Only time_period nor iteration_period can be specified for the OutputTimeControler {:s}".format(self))
        self.__time_period = time_period
        self.__iteration_period = iteration_period
        self.__next_output_time = self.__time_period
        self.__next_output_iteration = self.__iteration_period

    def __str__(self):
        return self.__class__.__name__ + " of the database {:s}".format(self.__id)

    def db_has_to_be_updated(self, time, iteration):
        """
        Return True if the iteration or time requires to write fields in the database
        and update the next output time or iteration
        """
        answer = False
        try:
            if time > self.__next_output_time:
                answer = True
                self.__next_output_time += self.__time_period
        except TypeError:
            if iteration > self.__next_output_iteration:
                answer = True
                self.__next_output_iteration += self.__iteration_period
        return answer


class OutputManager(object):
    """
    The manager of all outputs
    """
    __metaclass__ = Singleton

    def __init__(self):
        self.__db_build_infos = {}

    def register_database_time_ctrl(self, database_name, database_obj, delta_t):
        """
        Add a database to the manager. The database will be updated every deltat_t seconds
        """
        time_ctrl = OutputTimeControler(database_name, time_period=delta_t)
        dbinfos = DatabaseBuildInfos(database_obj, fields=[], time_controler=time_ctrl)
        self.__db_build_infos[database_name] = dbinfos
    
    def register_database_iteration_ctrl(self, database_name, database_obj, delta_it):
        """
        Add a database to the manager. The database will be updated every delta_it iterations
        """
        time_ctrl = OutputTimeControler(database_name, iteration_period=delta_it)
        dbinfos = DatabaseBuildInfos(database_obj, fields=[], time_controler=time_ctrl)
        self.__db_build_infos[database_name] = dbinfos

    def register_field(self, field_name, field_support, field_attr_name, *database_names):
        """
        Add a field to the manager. Each field will be stored in every database in arguments
        if specified else in every database registered
        """
        if database_names:
            for db in database_names:
                self.__db_build_infos[db].fields.append(
                        Field(name=field_name, owner=field_support, attr_name=field_attr_name))
        else:
            for build_infos in self.__db.values():
                build_infos.fields.append(Field(name=field_name, owner=field_support, attr_name=field_attr_name))

    def update(self, time, iteration):
        """
        If the current time given in argument is above the time of next output then
        the manager asks each of its database to save fields. The same appends for 
        a current iteration above the iteration of next output
        """
        for build_infos in self.__db_build_infos.values():
            if build_infos.time_controler.db_has_to_be_updated(time, iteration):
                build_infos.database_object.add_time(time)
                for field in build_infos.fields:
                    build_infos.database_object.add_field(field.name, getattr(field.owner, field.attr_name))

    def finalize(self):
        """
        Close all the database
        """
        for build_infos in self.__db_build_infos.values():
            build_infos.database_object.close()

    def __str__(self):
        msg = self.__class__.__name__
        msg += " manages following database: " + os.linesep
        for db_name in self.__db_build_infos:
            msg += "|_> " + db_name + os.linesep
            msg += "    Fields : {:s}".format(
                    ",".join([field.name for field in self.__db_build_infos[db_name].fields])) + os.linesep
        return msg
