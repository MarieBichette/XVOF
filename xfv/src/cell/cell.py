# -*- coding: utf-8 -*-
"""
Implementing class cell
"""

from abc import abstractmethod
from copy import deepcopy
import os
import numpy as np

from xfv.src.data.data_container import DataContainer
from xfv.src.fields.field import Field
from xfv.src.fields.fieldsmanager import FieldManager


class Cell:  # pylint: disable=too-many-public-methods, too-many-instance-attributes
    """
    A Cell object represents all the mesh cells.
    Its different members are, for most of them, numpy 1D-array of nbr_of_cells length.

    Memory layout is the same as in C/C++, i-e 'row wise'.
    """

    @classmethod
    def get_coordinates(cls, nbr_cells, topology, x_coord,
                        y_coord=None, z_coord=None):  # pylint: disable=too-many-arguments
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

        for ielem in range(nbr_cells):
            nodes_index = topology.get_nodes_belonging_to_cell(ielem)
            vec_coord[ielem][0] = x_coord[nodes_index].mean()
            if topology.dimension == 2:
                vec_coord[ielem][1] = y_coord[nodes_index].mean()
            if topology.dimension == 3:
                vec_coord[ielem][2] = z_coord[nodes_index].mean()
        return vec_coord

    def __init__(self, nbr_of_cells: int):
        """
        Constructor of the array of cells

        :param nbr_of_cells: number of cells
        """
        self.data = DataContainer()  # pylint: disable=no-value-for-parameter
        self._nbr_of_cells = nbr_of_cells
        self._dt = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._size_t = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._size_t_plus_dt = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._mass = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._fields_manager = FieldManager()

        self.cell_in_target = np.zeros(self.number_of_cells, dtype='bool')
        self.cell_in_projectile = np.zeros(self.number_of_cells, dtype='bool')

        # Default initialization for target material ---------------------
        # hydro :
        material_data = self.data.material_target.initial_values
        self._fields_manager["Density"] = Field(
            self._nbr_of_cells, material_data.rho_init, material_data.rho_init)
        self._fields_manager["Pressure"] = Field(
            self._nbr_of_cells, material_data.pression_init, material_data.pression_init)
        self._fields_manager["Energy"] = Field(
            self._nbr_of_cells, material_data.energie_init, material_data.energie_init)
        self._fields_manager["Pseudo"] = Field(self._nbr_of_cells)
        self._fields_manager["SoundVelocity"] = Field(self._nbr_of_cells)
        # Elasticity
        self._stress = np.zeros([self.number_of_cells, 3], dtype=np.float64, order='C')
        self._fields_manager["ShearModulus"] = Field(
            self._nbr_of_cells, material_data.shear_modulus_init, material_data.shear_modulus_init)
        self._fields_manager["YieldStress"] = Field(
            self._nbr_of_cells, material_data.yield_stress_init, material_data.yield_stress_init)
        # Porosity
        self._fields_manager["Porosity"] = Field(
            self._nbr_of_cells, material_data.porosity_init, material_data.porosity_init)

        # Coordonnées des mailles selon l'axe x
        self._coordinates_x = np.zeros(self.number_of_cells, dtype=float)

    def initialize_cell_fields(self, mask_node_target, mask_node_projectile, topology):
        """
        Initialisation of the cell fields and attributes of cell_in_target and cell_in_projectile

        :param mask_node_target: bool array for nodes in the target
        :param mask_node_projectile: bool array for nodes in the target
        :param topology: mesh connectivity object


        :type mask_node_target: numpy.array([nbr_of_cells, 1], dtype=bool, order='C')
        :type mask_node_projectile: numpy.array([nbr_of_cells, 1], dtype=bool, order='C')
        :type topology: Topology
        """
        # Part : mask_target
        node_indexes = np.where(mask_node_target)[0]
        cell_target = np.unique(
            topology.get_cells_in_contact_with_node(node_indexes)[1:-1].flatten())
        # [1:-1] => elimination of the extremal nodes because their connectivity is incomplete
        self.cell_in_target[cell_target] = True

        # Part : mask_projectile
        node_indexes = np.where(mask_node_projectile)[0]
        cell_projectile = np.unique(
            topology.get_cells_in_contact_with_node(node_indexes)[1:-1].flatten())
        # [1:-1] => elimination of the extremal nodes because their connectivity is incomplete
        self.cell_in_projectile[cell_projectile] = True

        try:
            print("Cells in the projectile : de {:} a {:}" .format(
                np.where(self.cell_in_projectile)[0][0], np.where(self.cell_in_projectile)[0][-1]))
            print("Cells in the target : de {:} a {:}".format(
                np.where(self.cell_in_target)[0][0], np.where(self.cell_in_target)[0][-1]))
        except IndexError:
            # case where no projectile or target exists.
            # self.cell_in_projectile or self.cell_in_target arrays are empty
            pass

        # Initialisation of the projectile material data if a projectile exists in XDATA.json
        # (mask is false everywhere if no projectile in the DataContainer
        # => transparent operation if no projectile declared in the DC)
        if self.data.data_contains_a_projectile:
            material_data = self.data.material_projectile.initial_values
            self.density.current_value[self.cell_in_projectile] = material_data.rho_init
            self.density.new_value[self.cell_in_projectile] = material_data.rho_init
            self.pressure.current_value[self.cell_in_projectile] = material_data.pression_init
            self.pressure.new_value[self.cell_in_projectile] = material_data.pression_init
            self.energy.current_value[self.cell_in_projectile] = material_data.energie_init
            self.energy.new_value[self.cell_in_projectile] = material_data.energie_init
            self.shear_modulus.current_value[self.cell_in_projectile] = \
                material_data.shear_modulus_init
            self.shear_modulus.new_value[self.cell_in_projectile] = \
                material_data.shear_modulus_init
            self.yield_stress.current_value[self.cell_in_projectile] = \
                material_data.yield_stress_init
            self.yield_stress.new_value[self.cell_in_projectile] = material_data.yield_stress_init
            self.porosity.current_value[self.cell_in_projectile] = material_data.porosity_init
            self.porosity.new_value[self.cell_in_projectile] = material_data.porosity_init

    @property
    def dt(self):  # pylint: disable=invalid-name
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
    def mass(self):
        """
        Mass of the cells
        """
        return self._mass

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
    def shear_modulus(self):
        """
        Elastic shear modulus
        """
        return self._fields_manager["ShearModulus"]

    @property
    def yield_stress(self):
        """
        Yield stress separating elastic from plastic behavior
        """
        return self._fields_manager["YieldStress"]

    @property
    def stress(self):
        """
        Cauchy stress tensor in the cells
        """
        return self._stress

    @property
    def stress_xx(self):
        """
        Cauchy stress tensor in the cells. 1D : component xx
        """
        return self._stress[:, 0]

    @property
    def stress_yy(self):
        """
        Stress tensor field : sigma_yy
        """
        return self._stress[:, 1]

    @property
    def stress_zz(self):
        """
        Stress tensor field : sigma_zz
        """
        return self._stress[:, 2]

    @property
    def porosity(self):
        """
        Porosity
        """
        return self._fields_manager["Porosity"]

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

    @property
    def coordinates_x(self):
        """
        Returns the coordinates of the center of the cell
        """
        return self._coordinates_x

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
        message += "==> sound velocity at t = {}".format(
            self.sound_velocity.current_value) + os.linesep
        message += "==> sound velocity at t+dt = {}".format(self.sound_velocity.new_value) + os.linesep
        message += "==> porosity at t = {}".format(self.porosity.current_value) + os.linesep
        message += "==> porosity at t+dt = {}".format(self.porosity.new_value) + os.linesep
        print(message)

    def increment_variables(self):
        """
        Variables incrementation
        """
        self._fields_manager.increment_fields()
        self._size_t[:] = self._size_t_plus_dt[:]

    @abstractmethod
    def compute_new_pressure(self, mask, delta_t):
        """
        Compute the pressure in the cells at time t + dt
        """

    @abstractmethod
    def compute_size(self, topology, node_coord):
        """
        Compute the size of the cells

        :param topology: topology of the mesh
        :param node_coord: array of nodal coordinates

        :type topology: Topology
        :type node_coord: numpy.array([nbr_of_nodes, 1], dtype=np.float64, order='C')
        """

    @abstractmethod
    def compute_new_size(self, *args, **kwargs):
        """
        Compute the new size of the cells
        """

    @abstractmethod
    def compute_mass(self):
        """
        Compute mass of the cells
        """

    @abstractmethod
    def compute_new_density(self, mask):
        """
        Compute the new density in the cells

        :param mask: boolean array to identify cells to be computed

        :type mask: np.array([nbr_of_cells, 1], dtype=bool)
        """

    @abstractmethod
    def compute_new_pseudo(self, time_step, mask):
        """
        Compute the new value of artificial viscosity in the cells

        :param time_step: time step
        :param mask: boolean array to identify cells to be computed

        :type time_step: float
        :type mask: np.array([nbr_of_cells, 1], dtype=bool)
        """

    @abstractmethod
    def compute_new_time_step(self, mask):
        """
        Compute the new value of critical time step in the cells

        :param mask: boolean array to identify cells to be computed

        :type mask: np.array([nbr_of_cells, 1], dtype=bool)
        """

    @abstractmethod
    def compute_new_porosity(self, time_step, porosity_model, mask):
        """
        Compute the new porosity according to the porosity model in XDATA
        :param time_step: float
        :param porosity_model: porosity model to compute 
        :type mask: np.array([nbr_of_cells, 1], dtype=bool)
        """

    @abstractmethod
    def compute_new_coordinates(self, topology, x_coord):
        """
        Compute the coordinates of the cell center
        :param topology: mesh nodal connectivity
        :param x_coord: coordinates of the nodes
        :param cell_size: size of the cells
        :return:
        """
