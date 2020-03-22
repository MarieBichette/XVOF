# -*- coding: utf-8 -*-
"""
Implements the TimeData, CellTimeData and NodeTimeData classes

:todo: find where the add_time_step_fields method should be called
:todo: self.__data should be a numpy array and thus there is no need for history_array property
:todo: in the add_time_step_fields use *args instead of a list of cell fields.
       Then remove the pylint pragma.
:todo: maybe used a ContextManager
"""
import os
from abc import ABCMeta, abstractmethod
import numpy as np


class TimeData(object):
    """
    This class is an abstract base class for CellTimeData and NodeTimeData classes.
    """
    __metaclass__ = ABCMeta

    def __init__(self, attribute, id_number, o_path):
        """
        :param id_number: item id
        :type id_nubmer: int
        :param o_path: output file path
        :type o_path: str
        """
        self._item_id = id_number
        self._path = "./{:s}/{:s}_history_{:<03d}.dat".format(o_path, attribute, id_number)
        try:
            self._fd = open(self._path, "w")
        except IOError:
            raise SystemExit("Unable to open the file under : {:s}".format(self._path))

    @abstractmethod
    def add_time_step_fields(self, *args):
        """
        Stocks the fields from the current time step
        """
        pass

    @abstractmethod
    def write_fields_history(self):
        """
        Write the history data in appropriate file
        """
        pass

    def close_file(self):
        """
        Close the output file
        """
        self._fd.close()


class CellTimeData(TimeData):
    """
    This class writes a file which holds the time evolution of the fields of a cell
    """
    def __init__(self, id_number, o_path):
        """
        :param id_number: cell id
        :param o_path: output file path
        """
        super(CellTimeData, self).__init__("cell", id_number, o_path)
        self._data = []

    @property
    def header(self):
        """
        Returns the header for the cell file

        :return: the header for the cell file
        """
        msg = "---Time History for cell :" + str(self._item_id) + os.linesep
        msg += "{:^15s}    {:^15s}    {:^15s}    {:^15s}    {:^15s}" .format(
            "Time[s]", "Density[kg/m3]", "Pressure[MPa]", "Internal Energy[J/m3]",
            "Artificial Viscosity[MPa]")
        msg += os.linesep
        return msg

    def add_time_step_fields(self, time, density, pressure, energy, pseudo):  # pylint:disable=too-many-arguments, arguments-differ
        """
        Append the history file with the current cell data

        :param time: float, real simulation time
        :param density: density field
        :param pressure: pressure field
        :param energy: internal energy field
        :param pseudo: artificial viscosity field

        :warning: this method seems to be unused
        """
        self._data.append(
            [time, density[self._item_id], pressure[self._item_id],
             energy[self._item_id], pseudo[self._item_id]])

    @property
    def history_array(self):
        """
        Returns the underlying array
        """
        return np.array(self._data)

    def write_fields_history(self):
        """
        Write the history data in appropriate file
        """
        np.savetxt(self._path, self.history_array,
                   fmt=['%+10.9e', '%+10.9e', '%+10.9e', '%+10.9e', '%+10.9e'],
                   header=self.header)


class NodeTimeData(TimeData):
    """
    This class writes a file which holds the time evolution of the fields of a node
    """
    def __init__(self, id_number, o_path):
        """
        :param id_number: node id
        :param o_path: output file path
        """
        super(NodeTimeData, self).__init__("node", id_number, o_path)
        self._data = []

    @property
    def header(self):
        """
        Returns the header for the node file

        :return: the header for the node file
        """
        msg = "---Time History for node :" + str(self._item_id) + os.linesep
        msg += "{:^15s}    {:^15s}    {:^15s}".format("Time[s]", "Position[m]", "Vitesse[m/s]")
        msg += os.linesep
        return msg

    def add_time_step_fields(self, time, position, velocity):  # pylint: disable=arguments-differ
        """
        Append the history file with the current node data

        :param time: float, real simulation time
        :param position: node position
        :param velocity: node velocity

        :warning: this method seems unused
        """
        self._data.append([time, position[self._item_id][0], velocity[self._item_id][0]])

    @property
    def history_array(self):
        """
        Returns the underlying array
        """
        return np.array(self._data)

    def write_fields_history(self):
        """
        Write the history data in appropriate file
        """
        np.savetxt(self._path, self.history_array,
                   fmt=['%+10.9e', '%+10.9e', '%+10.9e'],
                   header=self.header)
