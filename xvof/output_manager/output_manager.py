"""
Implementing the OutputManager class
"""
from collections import namedtuple
import os
from xvof.utilities.singleton import Singleton
from xvof.output_manager.outputtimecontroler import OutputTimeControler

DatabaseBuildInfos = namedtuple("DatabaseBuildInfos", ["database_object", "fields", "time_controler"])
Field = namedtuple("Field", ["name", "owner", "attr_name", "indices"])


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

    def register_field(self, field_name, field_support, field_attr_name, indices=None, database_names=None):
        """
        Add a field to the manager. Each field will be stored in every database in arguments
        if specified else in every database registered
        """
        if database_names is not None:
            for db in database_names:
                self.__db_build_infos[db].fields.append(
                    Field(name=field_name, owner=field_support, attr_name=field_attr_name, indices=indices))
        else:
            for build_infos in self.__db_build_infos.values():
                build_infos.fields.append(Field(name=field_name, owner=field_support,
                                                attr_name=field_attr_name, indices=indices))

    def update(self, time, iteration):
        """
        If the current time given in argument is above the time of next output then
        the manager asks each of its database to save fields. It's the same for
        a current iteration above the iteration of next output
        """
        for build_infos in self.__db_build_infos.values():
            if build_infos.time_controler.db_has_to_be_updated(time, iteration):
                build_infos.database_object.add_time(time)
                for field in build_infos.fields:
                    if field.indices is None:
                        build_infos.database_object.add_field(field.name, getattr(field.owner, field.attr_name))
                    else:
                        build_infos.database_object.add_field(
                            field.name, getattr(field.owner, field.attr_name).__getitem__(field.indices))


    def finalize(self):
        """
        Close all the database
        """
        print("Flushing and writing database outputs!")
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
