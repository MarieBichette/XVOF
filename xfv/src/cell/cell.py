# -*- coding: iso-8859-1 -*-
"""
Implementing class cell
"""

import numpy as np
from abc import abstractmethod
from copy import deepcopy
import os

from xfv.src.data.data_container import DataContainer
from xfv.src.fields.field import Field
from xfv.src.fields.fieldsmanager import FieldManager


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
            nodes_index = topology.get_nodes_belonging_to_cell(ielem)
            vec_coord[ielem][0] = x_coord[nodes_index].mean()
            if topology.dimension == 2:
                vec_coord[ielem][1] = y_coord[nodes_index].mean()
            if topology.dimension == 3:
                vec_coord[ielem][2] = z_coord[nodes_index].mean()
        return vec_coord

    def __init__(self, nbr_of_cells):
        self._nbr_of_cells = nbr_of_cells
        self._dt = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._size_t = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._size_t_plus_dt = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._mass = np.zeros(self._nbr_of_cells, dtype=np.float64, order='C')
        self._fields_manager = FieldManager()

        self.cell_in_target = np.zeros(self.number_of_cells, dtype='bool')
        self.cell_in_projectile = np.zeros(self.number_of_cells, dtype='bool')

        # initialisation par d�faut avec material_target ---------------------
        # hydro :
        material_data = DataContainer().material_target.initial_values
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

    def initialize_cell_fields(self, mask_node_target, mask_node_projectile, topology):
        """
        Initialisation des champs aux mailles et des caract�ristiques cell_in_target et
        cell_in_projectile
        :param mask_node_target: tableau de bool pour les noeuds dans la cible
        :param mask_node_projectile: tableau de bool pour les noeuds dans le projectile
        :param topology: Topologie donnant les tableaux de connectivit�
        :return:
        """
        # Partie mask_cible
        indice_noeuds = np.where(mask_node_target)[0]
        cell_target = np.unique(topology.get_cells_in_contact_with_node(indice_noeuds)[1:-1].flatten())
        # on �limine les noeuds extr�mes :
        # 1 pour ne pas prendre la derni�re maille du projectile qui est connect�e
        # aussi � un noeud target -1 pour ne pas prendre la maille connect�ee au dernier
        # noeud de la cible, qui n'existe pas
        self.cell_in_target[cell_target] = True

        # Partie mask_projectile
        indice_noeuds = np.where(mask_node_projectile)[0]
        cell_projectile = np.unique(topology.get_cells_in_contact_with_node(indice_noeuds)[1:-1].flatten())
        # on �limine les noeuds extr�mes :
        # -1 pour ne pas prendre la premi�re maille de la cible qui est aussi
        # connect�e � un noeud projectile 1 pour ne pas prendre la maille connect�ee au
        # premier noeud du projectile, qui n'existe pas
        self.cell_in_projectile[cell_projectile] = True

        try:
            print("Cells in the projectile : de {:} a {:}" .format(
                np.where(self.cell_in_projectile)[0][0], np.where(self.cell_in_projectile)[0][-1]))
            print("Cells in the target : de {:} a {:}".format(
                np.where(self.cell_in_target)[0][0], np.where(self.cell_in_target)[0][-1]))
        except IndexError:
            # pour g�rer  les exceptions o� il n'y a pas de projectile ou de target
            # (cas tableau vide [index])
            pass

        # correction de l'init si materiau projectile d�clar� dans XDATA.xml
        # (le mask est vide si pas de projectile donc transparent quand il n'y a pas
        # de projectile d�clar�)
        material_data = DataContainer().material_projectile.initial_values
        self.density.current_value[self.cell_in_projectile] = material_data.rho_init
        self.density.new_value[self.cell_in_projectile] = material_data.rho_init
        self.pressure.current_value[self.cell_in_projectile] = material_data.pression_init
        self.pressure.new_value[self.cell_in_projectile] = material_data.pression_init
        self.energy.current_value[self.cell_in_projectile] = material_data.energie_init
        self.energy.new_value[self.cell_in_projectile] = material_data.energie_init
        self.shear_modulus.current_value[self.cell_in_projectile] = material_data.shear_modulus_init
        self.shear_modulus.new_value[self.cell_in_projectile] = material_data.shear_modulus_init
        self.yield_stress.current_value[self.cell_in_projectile] = material_data.yield_stress_init
        self.yield_stress.new_value[self.cell_in_projectile] = material_data.yield_stress_init


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
        message += "==> sound velocity at t = {}".format(
            self.sound_velocity.current_value) + os.linesep
        message += "==> sound velocity at t+dt = {}".format(self.sound_velocity.new_value)
        print(message)

    def increment_variables(self):
        """
        Variables incrementation
        """
        self._fields_manager.incrementFields()
        self._size_t[:] = self._size_t_plus_dt[:]

    @abstractmethod
    def compute_new_pressure(self, mask, dt):
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
    def compute_mass(self):
        """
        Compute mass of the cells
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
