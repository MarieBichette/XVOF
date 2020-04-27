# -*- coding: iso-8859-1 -*-
"""
Classe d�finissant le gestionnaire d'images

@todo: H�riter de la classe Figure
@todo: translate in english
@todo: incorporate the calculation of the time at which images have to be shown (move from Vnr1D...py)
"""
import matplotlib.pyplot as plt
import numpy as np
import os
from collections import namedtuple
from .physic_figure import PhysicFigure
from xfv.src.utilities.singleton import Singleton
from xfv.src.output_manager.outputtimecontroler import OutputTimeControler

max_coord = 0.015
Field = namedtuple("Field", ["label", "titre", "val_min", "val_max", "results_path"])

PressureField = Field("Pression [Pa]", "Champ de pression", -20e+09, 20e+09,
                      "./RESULTATS/PressureField")
DensityField = Field("Masse volumique [kg/m3]", "Champ de densite", 7000.0, 12000.0,
                     "./RESULTATS/DensityField")
InternalEnergyField = Field("Energie interne [J/kg]", "Champ d energie interne", 0, 100000.0,
                            "./RESULTATS/InternalEnergyField")
PseudoViscosityField = Field("Pseudoviscosite [Pa]", "Champ de pseudoviscosite", 0, 1.5e+09,
                             "./RESULTATS/PseudoViscosityField")
CellPositionField = Field("Position [m]", "Champ de position", 0.0, max_coord,
                          "./RESULTATS/CellPositionField")
NodePositionField = Field("Position [m]", "Champ de position", 0.0, max_coord,
                          "./RESULTATS/NodePositionField")
NodeVelocityField = Field("Vitesse [m/s]", "Champ de vitesse", -1000.0, 1000.0,
                          "./RESULTATS/NodeVelocityField")
DevStressField = Field("Stress Dev [Pa]", "Champ de dev stress", -1.5e+09,
                        1.5e+09, "./RESULTATS/DevStressField")
StressXXField = Field("Stress XX [Pa]", "Champ de stress XX", -20e+09, 20e+09,
                      "./RESULTATS/StressXXField")


class FigureManager(object, metaclass=Singleton):
    """
    Gestionnaire de figures
    """

    def __init__(self, mesh_instance, dump=False):
        self.__mesh_instance = mesh_instance
        self.__figures_mailles = []
        self.__figures_noeuds = []
        self.__champs_mailles = {}
        self.__champs_noeuds = {}
        self.update_fields()
        self.__time_ctrl = None
        self.__dump = dump
        self.interface = int(np.where(self.__mesh_instance.cells.cell_in_target)[0][0])

    def set_time_controler(self, deltat_t):
        """
        The figures will be updated every delta_t seconds

        :param deltat_t: interval between two figures
        """
        self.__time_ctrl = OutputTimeControler(identifier="FigureManager", time_period=deltat_t)

    def set_iteration_controler(self, deltat_it):
        """
        The figures will be updated every delta_it iterations

        :param deltat_it: interval between two figures
        """
        self.__time_ctrl = OutputTimeControler(identifier="FigureManager",
                                               iteration_period=deltat_it)

    def update_fields(self):
        """ MAJ des champs par appel des propri�t�s du maillage"""
        self.__champs_mailles = \
             {CellPositionField: self.__mesh_instance.cells_coordinates[:],
             PressureField: self.__mesh_instance.pressure_field,
             # DensityField: self.__mesh_instance.density_field,
             # DevStressField: self.__mesh_instance.deviatoric_stress_field,
             StressXXField: self.__mesh_instance.stress_xx_field,
             }
        self.__champs_noeuds = \
            {
             NodePositionField: self.__mesh_instance.nodes_coordinates[:],
             NodeVelocityField: self.__mesh_instance.velocity_field
             }

    def create_figure_for_cell_field(self, field_X, field_Y):
        """
        Cr�ation des figures pour les champs aux mailles
        (l'axe des X est donc l'abscisse des mailles)
        """
        try:
            X = self.__champs_mailles[field_X]
        except ValueError as ve:
            print("Le champ {} est inconnu!".format(field_X))
            raise ve
        try:
            Y = self.__champs_mailles[field_Y]
        except ValueError as ve:
            print("Le champ {} est inconnu!".format(field_Y))
            raise ve
        if self.__dump:
            phyfig = PhysicFigure(X, Y, xlabel=field_X.label, ylabel=field_Y.label,
                                  titre=field_Y.titre, interface_id=self.interface,
                                  save_path=field_Y.results_path)
        else:
            phyfig = PhysicFigure(X, Y, xlabel=field_X.label, ylabel=field_Y.label,
                                  titre=field_Y.titre, interface_id=self.interface)
        phyfig.set_y_limit(field_Y.val_min, field_Y.val_max)
        phyfig.set_x_limit(field_X.val_min, field_X.val_max)
        return phyfig

    def create_figure_for_node_field(self, field_X, field_Y):
        """
        Cr�ation des figures pour les champs aux noeuds
        (l'axe des X est donc l'abscisse des noeuds)
        """
        try:
            X = np.array(self.__champs_noeuds[field_X])
        except ValueError as ve:
            print("Le champ {} est inconnu!".format(field_X))
            raise ve
        try:
            Y = np.array(self.__champs_noeuds[field_Y])
        except ValueError as ve:
            print("Le champ {} est inconnu!".format(field_Y))
            raise ve
        if self.__dump:
            phyfig = PhysicFigure(X, Y, xlabel=field_X.label, ylabel=field_Y.label,
                                  titre=field_Y.titre, interface_id=self.interface + 1,
                                  save_path=field_Y.results_path)
        else:
            phyfig = PhysicFigure(X, Y, xlabel=field_X.label, ylabel=field_Y.label,
                                  titre=field_Y.titre, interface_id=self.interface + 1)
        phyfig.set_y_limit(field_Y.val_min, field_Y.val_max)
        phyfig.set_x_limit(field_X.val_min, field_X.val_max)
        return phyfig

    def populate_figs(self):
        """
        Cr�ation des figures associ�es � chacun des champs et ajout �
        la liste des figures
        """
        champ_X = CellPositionField
        for champ_Y in list(self.__champs_mailles.keys()):
            if champ_Y != champ_X:
                fig = self.create_figure_for_cell_field(champ_X, champ_Y)
                self.__figures_mailles.append((fig, champ_X, champ_Y))
        champ_X = NodePositionField
        for champ_Y in list(self.__champs_noeuds.keys()):
            if champ_Y != champ_X:
                fig = self.create_figure_for_node_field(champ_X, champ_Y)
                self.__figures_noeuds.append((fig, champ_X, champ_Y))
        plt.show(block=False)

        if self.__dump:
            self.create_reps()

    def update_figs(self, title_compl=None):
        """
        MAJ des champs puis des figures correspondantes
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
                msg+= "\n Abandon to avoid to override data"
                raise SystemExit(msg)
            else:
                os.makedirs(path)