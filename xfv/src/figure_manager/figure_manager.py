# -*- coding: utf-8 -*-
# pylint: disable=invalid-name

"""
A class to manage the figures and animations

@todo: incorporate the calculation of the time at which images have to be shown
(move from Vnr1D...py)
"""
from collections import namedtuple
import os
import matplotlib.pyplot as plt
import numpy as np
from xfv.src.figure_manager.physic_figure import PhysicFigure
from xfv.src.utilities.singleton import Singleton
from xfv.src.output_manager.outputtimecontroler import OutputTimeControler

max_coord = 0.015
Field = namedtuple("Field", ["label", "titre", "val_min", "val_max", "results_path"])

PressureField = Field("Pressure [Pa]", "Pressure Field", -20e+09, 20e+09, "./RESULTS/PressureField")
DensityField = Field("Density [kg/m3]", "Density Field", 7000.0, 12000.0, "./RESULTS/DensityField")
InternalEnergyField = Field("Internal Energy [J/kg]", "Internal Energy field", 0, 100000.0,
                            "./RESULTS/InternalEnergyField")
PseudoViscosityField = Field("Artificial viscosity [Pa]", "Artificial viscosity field", 0, 1.5e+09,
                             "./RESULTS/PseudoViscosityField")
CellPositionField = Field("Position [m]", "Cell position field", 0.0, max_coord,
                          "./RESULTS/CellPositionField")
NodePositionField = Field("Position [m]", "Node position field", 0.0, max_coord,
                          "./RESULTS/NodePositionField")
NodeVelocityField = Field("Velocity [m/s]", "Velocity field", -1000.0, 1000.0,
                          "./RESULTS/NodeVelocityField")
DevStressField = Field("Stress Dev [Pa]", "Deviatoric Stress XX field", -1.5e+09, 1.5e+09,
                       "./RESULTS/DevStressField")
StressXXField = Field("Stress XX [Pa]", "Stress XX field", -2e+10, 2e+10, "./RESULTS/StressXXField")


class FigureManager(metaclass=Singleton):
    """
    A manager for figures and animations during calculation
    """

    def __init__(self, mesh_instance, dump=False):
        """
        Creation of the instance

        :param mesh_instance:
        :param dump:
        """
        self.__mesh_instance = mesh_instance
        self.__figures_mailles = []
        self.__figures_noeuds = []
        self.__champs_mailles = {}
        self.__champs_noeuds = {}
        self.update_fields()
        self.__time_ctrl = None
        self.__dump = dump
        self.interface = int(np.where(self.__mesh_instance.cells.cell_in_target)[0][0])

    def set_time_controler(self, deltat_t: float):
        """
        The figures will be updated every delta_t seconds

        :param deltat_t: interval between two figures
        """
        self.__time_ctrl = OutputTimeControler(identifier="FigureManager", time_period=deltat_t)

    def set_iteration_controler(self, deltat_it: float):
        """
        The figures will be updated every delta_it iterations

        :param deltat_it: interval between two figures
        """
        self.__time_ctrl = OutputTimeControler(identifier="FigureManager",
                                               iteration_period=deltat_it)

    def update_fields(self):
        """
        Update the fields to plot using the mesh update methods
        """
        # Add / remove here the fields to be plotted
        self.__champs_mailles = {CellPositionField: self.__mesh_instance.cells_coordinates[:],
                                 PressureField: self.__mesh_instance.pressure_field,
                                 StressXXField: self.__mesh_instance.stress_xx_field
                                 }
        self.__champs_noeuds = {NodePositionField: self.__mesh_instance.nodes_coordinates[:],
                                NodeVelocityField: self.__mesh_instance.velocity_field
                                }

    def create_figure_for_cell_field(self, field_x, field_y):
        """
        Creation of figures for cell fields. Abscissa is thus the cell coordinates

        :param field_x: abscissa of the figure
        :param field_y: ordinate of the figure
        """
        try:
            x_value = self.__champs_mailles[field_x]
        except ValueError as ve:
            print("Field {} is unknown !".format(field_x))
            raise ve
        try:
            y_value = self.__champs_mailles[field_y]
        except ValueError as ve:
            print("Field {} is unknown !".format(field_y))
            raise ve
        if self.__dump:
            phyfig = PhysicFigure(x_value, y_value, xlabel=field_x.label, ylabel=field_y.label,
                                  titre=field_y.titre, interface_id=self.interface,
                                  save_path=field_y.results_path)
        else:
            phyfig = PhysicFigure(x_value, y_value, xlabel=field_x.label, ylabel=field_y.label,
                                  titre=field_y.titre, interface_id=self.interface)
        phyfig.set_y_limit(field_y.val_min, field_y.val_max)
        phyfig.set_x_limit(field_x.val_min, field_x.val_max)
        return phyfig

    def create_figure_for_node_field(self, field_x, field_y):
        """
        Creation of figures for nodal fields. Abscissa is thus the node coordinates

        :param field_x: abscissa of the figure
        :param field_y: ordinate of the figure
        """
        try:
            x_value = np.array(self.__champs_noeuds[field_x])
        except ValueError as ve:
            print("Le champ {} est inconnu!".format(field_x))
            raise ve
        try:
            y_value = np.array(self.__champs_noeuds[field_y])
        except ValueError as ve:
            print("Le champ {} est inconnu!".format(field_y))
            raise ve
        if self.__dump:
            phyfig = PhysicFigure(x_value, y_value, xlabel=field_x.label, ylabel=field_y.label,
                                  titre=field_y.titre, interface_id=self.interface + 1,
                                  save_path=field_y.results_path)
        else:
            phyfig = PhysicFigure(x_value, y_value, xlabel=field_x.label, ylabel=field_y.label,
                                  titre=field_y.titre, interface_id=self.interface + 1)
        phyfig.set_y_limit(field_y.val_min, field_y.val_max)
        phyfig.set_x_limit(field_x.val_min, field_x.val_max)
        return phyfig

    def populate_figs(self):
        """
        Creation of the figures associtated with fields and add to the list of figures
        """
        champ_x = CellPositionField
        for champ_y in list(self.__champs_mailles.keys()):
            if champ_y != champ_x:
                fig = self.create_figure_for_cell_field(champ_x, champ_y)
                self.__figures_mailles.append((fig, champ_x, champ_y))
        champ_x = NodePositionField
        for champ_y in list(self.__champs_noeuds.keys()):
            if champ_y != champ_x:
                fig = self.create_figure_for_node_field(champ_x, champ_y)
                self.__figures_noeuds.append((fig, champ_x, champ_y))
        plt.show(block=False)

        if self.__dump:
            self.create_reps()

    def update_figs(self, title_compl=None):
        """
        Update the fields and the associated figures

        :param title_compl: title of the figure
        """
        self.update_fields()
        for (fig, champ_x, champ_y) in self.__figures_mailles:
            champ_x_valeurs = self.__champs_mailles[champ_x]
            champ_y_valeurs = self.__champs_mailles[champ_y]
            fig.update(champ_x_valeurs, champ_y_valeurs, title_compl)

        for (fig, champ_x, champ_y) in self.__figures_noeuds:
            champ_x_valeurs = self.__champs_noeuds[champ_x]
            champ_y_valeurs = self.__champs_noeuds[champ_y]
            fig.update(champ_x_valeurs, champ_y_valeurs, title_compl)

        plt.show(block=False)

    def update(self, time, iteration):
        """
        If the current time given in argument is above the time of next output then
        the manager asks each of its database to save fields. It's the same for
        a current iteration above the iteration of next output

        :param time: simulation time
        :param iteration: id of the current iteration
        """
        if self.__time_ctrl is not None:
            if self.__time_ctrl.db_has_to_be_updated(time, iteration):
                self.update_figs("t={:5.4g} us".format(time / 1.e-06))

    def create_reps(self):
        """
        Creation of the reps where data is saved
        """
        for (_, _, field) in self.__figures_mailles + self.__figures_noeuds:
            path = field.results_path
            if os.path.exists(path):
                msg = "Path {:} already exists !".format(path)
                msg += "\n Abandon to avoid to override data"
                raise SystemExit(msg)
            # else
            os.makedirs(path)
