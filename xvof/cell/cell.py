# -*- coding: iso-8859-1 -*-
"""
Implementing class cell
"""

import numpy as np
from abc import abstractmethod
from copy import deepcopy
import os

from xvof.data.data_container import DataContainer
from xvof.fields.field import Field
from xvof.fields.fieldsmanager import FieldManager


class Cell(object):
    """
    A Cell object represents all the mesh cells.
    Its different members are, for most of them, numpy 1D-array of nbr_of_cells length.

    Memory layout is the same as in C/C++, i-e 'row wise'. 
    """

    @classmethod
    def get_coordinates(cls, nbr_cells, topology, x_coord, y_coord=None,  z_coord=None):
        """
        Return the vector of cell center coordinates at time t

        :param nbr_cells: number of cells
        :param topology: topology 
        :param x_coord: x coordinate vector
        :param y_coord: y coordinate vector
        :param z_coord: z coordinate vector

        :type topology: Topology
        :type x_coord: numpy.array([nbr_of_cells, 1], dtype=np.float64, order='C')
        :type y_coord: numpy.array([nbr_of_cells, 1], dtype=np.float64, order='C')
        :type z_coord: numpy.array([nbr_of_cells, 1], dtype=np.float64, order='C')
        :return: the vector of cell center coordinates at time t
        :rtype: numpy.array([nbr_of_cells, topology.dimension], dtype=np.float64, order='C')
        """
        vec_coord = np.zeros([nbr_cells, topology.dimension])

        for ielem in xrange(nbr_cells):
            nodes_index = topology.getNodesBelongingToCell(ielem)
            vec_coord[ielem][0] = x_coord[nodes_index].mean()
            if topology.dimension == 2:
                vec_coord[ielem][1] = y_coord[nodes_index].mean()
            if topology.dimension == 3:
                vec_coord[ielem][2] =  z_coord[nodes_index].mean()
        return vec_coord

    def __init__(self, nbr_of_cells):
        self._nbr_of_cells = nbr_of_cells
        self._dt = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._size_t = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._size_t_plus_dt = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._fields_manager = FieldManager()
        self._fields_manager["Density"] = Field(self._nbr_of_cells, DataContainer().material.rho_init,
                                                DataContainer().material.rho_init)
        self._fields_manager["Pressure"] = Field(self._nbr_of_cells, DataContainer().material.pression_init,
                                                 DataContainer().material.pression_init)
        self._fields_manager["Energy"] = Field(self._nbr_of_cells, DataContainer().material.energie_init,
                                               DataContainer().material.energie_init)
        self._fields_manager["Pseudo"] = Field(self._nbr_of_cells)
        self._fields_manager["SoundVelocity"] = Field(self._nbr_of_cells)

    @property
    def dt(self):
        """
        Critical time step in cells
        """
        return self._dt

    @property
    def size_t(self):
        """
        Size (length, area, volume) of the cells at time t
        """
        return self._size_t

    @property
    def size_t_plus_dt(self):
        """
        Size (length, area, volume) of the cells at time t + dt
        """
        return self._size_t_plus_dt

    @property
    def density(self):
        """
        Density in the cells
        """
        return self._fields_manager['Density']

    @property
    def pressure(self):
        """
        Pressure in the cells
        """
        return self._fields_manager['Pressure']

    @property
    def sound_velocity(self):
        """
        Sound velocity in the cells
        """
        return self._fields_manager['SoundVelocity']

    @property
    def energy(self):
        """
        Internal energy in the cells
        """
        return self._fields_manager['Energy']

    @property
    def pseudo(self):
        """
        Artificial viscosity in the cells
        """
        return self._fields_manager['Pseudo']

    @property
    def fields_manager(self):
        """
        Return a copy of the field manager
        """
        return deepcopy(self._fields_manager)

    @property
    def number_of_cells(self):
        """
        Number of cells
        """
        return self._nbr_of_cells

    def __str__(self):
        message = "Number of cells: {:d}".format(self.number_of_cells)
        return message

    def print_infos(self):
        """
        Print the fields in the cells
        """
        message = os.linesep + "{:s} ".format(self.__class__.__name__) + os.linesep
        message += "==> number of cells = {:d}".format(self.number_of_cells) + os.linesep
        message += "==> size at t = {}".format(self.size_t) + os.linesep
        message += "==> size at t+dt = {}".format(self.size_t_plus_dt) + os.linesep
        message += "==> density at t = {}".format(self.density.current_value) + os.linesep
        message += "==> density at t+dt = {}".format(self.density.new_value) + os.linesep
        message += "==> pressure at t = {}".format(self.pressure.current_value) + os.linesep
        message += "==> pressure at t+dt = {}".format(self.pressure.new_value) + os.linesep
        message += "==> internal energy at t = {}".format(self.energy.current_value) + os.linesep
        message += "==> internal energy at t+dt = {}".format(self.energy.new_value) + os.linesep
        message += "==> sound velocity at t = {}".format(self.sound_velocity.current_value) + os.linesep
        message += "==> sound velocity at t+dt = {}".format(self.sound_velocity.new_value)
        print message

    def increment_variables(self):
        """
        Variables incrementation
        """
        self._fields_manager.incrementFields()
        self._size_t[:] = self._size_t_plus_dt[:]

    @abstractmethod
    def compute_new_pressure(self, mask):
        """
        Compute the pressure in the cells at time t + dt
        """

    @abstractmethod
    def compute_size(self, topologie, vecteur_coord_noeuds):
        """
        Compute the size of the cells
        """

    @abstractmethod
    def compute_new_size(self, *args, **kwargs):
        """
        Compute the new size of the cells
        """

    @abstractmethod
    def compute_new_density(self, mask):
        """
        Compute the new density in the cells
        """

    @abstractmethod
    def compute_new_pseudo(self, time_step, mask):
        """
        Compute the new value of artificial viscosity in the cells
        """

    @abstractmethod
    def compute_new_time_step(self, mask):
        """
        Compute the new value of critical time step in the cells
        """
