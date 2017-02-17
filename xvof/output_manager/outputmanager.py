"""
Implementing the OutputManager class
"""
from collections import namedtuple
import os
from xvof.utilities.singleton import Singleton
from xvof.output_manager.outputtimecontroler import OutputTimeControler

DatabaseBuildInfos = namedtuple("DatabaseBuildInfos", ["database_object", "fields", "time_controler"])
Field = namedtuple("Field", ["name", "owner", "attr_name", "indexes"])


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

    def register_field(self, field_name, field_support, field_attr_name, indexes=None, database_names=None):
        """
        Add a field to the manager. Each field will be stored in every database in arguments
        if specified else in every database registered
        """
        if database_names is not None:
            for db in database_names:
                self.__db_build_infos[db].fields.append(
                    Field(name=field_name, owner=field_support, attr_name=field_attr_name, indexes=indexes))
        else:
            for build_infos in self.__db_build_infos.values():
                build_infos.fields.append(Field(name=field_name, owner=field_support,
                                                attr_name=field_attr_name, indexes=indexes))

    def register_all_fields(self, cells, nodes, database_id, cell_indexes=None, node_indexes=None):
        """
        Add all fields to the manager.

        :param cells: cells from which fields must be printed
        :param nodes: nodes from which fields must be printed
        :param database_id: identifier of the database
        :param cell_indexes: indexes of the cells to be printed
        :param node_indexes: indexes of the nodes to be printed
        """
        if cell_indexes is None and node_indexes is None:
            cell_indexes = slice(0, cells.number_of_cells)
            node_indexes = slice(0, nodes.number_of_nodes)
        self.register_field("NodeStatus", nodes, ("enriched",), database_names=[database_id], indexes=node_indexes)
        self.register_field("NodeCoordinates", nodes, ("xt",), database_names=[database_id], indexes=node_indexes)
        self.register_field("ClassicalNodeVelocity", nodes, ("umundemi",), database_names=[database_id], indexes=node_indexes)
        self.register_field("EnrichedNodeVelocity", nodes, ("umundemi_enriched",), database_names=[database_id], indexes=node_indexes)
        self.register_field("CellStatus", cells, ("enriched",), database_names=[database_id], indexes=cell_indexes)
        self.register_field("CellSize", cells, ("size_t_plus_dt",), database_names=[database_id], indexes=cell_indexes)
        self.register_field("CellLeftSize", cells, ("left_size", "new_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("CellRightSize", cells, ("right_size", "new_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalPressure", cells, ("pressure", "new_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("EnrichedPressure", cells, ("pressure", "new_enr_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalDensity", cells, ("density", "new_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("EnrichedDensity", cells, ("density", "new_enr_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalInternalEnergy", cells, ("energy", "new_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("EnrichedInternalEnergy", cells, ("energy", "new_enr_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalSoundVelocity", cells, ("sound_velocity", "new_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("EnrichedSoundVelocity", cells, ("sound_velocity", "new_enr_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalArtificialViscosity", cells, ("pseudo", "new_value"), database_names=[database_id], indexes=cell_indexes)
        self.register_field("EnrichedArtificialViscosity", cells, ("pseudo", "new_enr_value"), database_names=[database_id], indexes=cell_indexes)

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
                    if field.indexes is not None:
                        value = getattr(field.owner, field.attr_name[0])
                        for attr_name in field.attr_name[1:]:
                            value = getattr(value, attr_name)
                        build_infos.database_object.add_field(
                            field.name, value.__getitem__(field.indexes), support=field.owner.__class__.__name__)

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
