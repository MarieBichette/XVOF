#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une fichier contenant les champs d'un noeud ou d'un élément en fonction du temps
Extraction des données d'un élément ou un noeud et enregistrement dans un fichier .dat
"""
import os
from abc import ABCMeta, abstractmethod
import numpy as np

class TimeData(object):
    """
    Une classe abstraite pour les objets de type TimeData
    """
    __metaclass__ = ABCMeta

    def __init__(self, attribute, id_number , o_path):
        """

        :param attribute:
        :param id_number:
        :param o_path:
        """
        self.__attribute = attribute
        self._item_id = id_number
        self._path = ("./{:s}/{:s}_history_{:<03d}.dat").format(o_path, attribute, id_number)
        try:
            self._fd = open(self._path,"w")
        except IOError:
            raise SystemExit("Unable to open the file under : {:s}".format(self._path))

    
    @abstractmethod
    def add_time_step_fields(self):
        """
        stocke les données du pas de temps actuel dans des listes
        """
        pass

    @abstractmethod
    def write_fields_history(self):
        """
        write the history data in appropriate file
        """
        pass

    def close_file(self):
        """
        Ferme le fichier associé
        """
        self._fd.close()


class CellTimeData(TimeData):
    """
    Classe définissant un jeu de donnée time history pour un élément
    """
    def __init__(self, id_number, o_path):
        """
        initialisation comme particularisation de TimeData
        :param id_number: numéro de l'élément
        :param o_path: chemin pour enregistrer dans le bon dossier
        """
        super(CellTimeData, self).__init__("cell", id_number, o_path)
        self._data = []


    @property
    def header(self):
        """
        crée le fichier pour cell time history
        :return:
        """
        msg = "---Time History for cell :" + str(self._item_id) + os.linesep
        msg += "{:^15s}    {:^15s}    {:^15s}    {:^15s}    {:^15s}" .format("Time[s]", "Density[kg/m3]", "Pressure[MPa]",
                                                                             "Internal Energy[J/m3]", "Artificial Viscosity[MPa]")
        msg += os.linesep
        return msg

    def add_time_step_fields(self, time, density, pressure, energy, pseudo):
        """
        Append the history file with the current cell data
        :var
        - path : str, path to the history data file
        - time : float, real simulation time
        - density, pressure, energy, pseudo : array, cell fields to be added
        """
        # self._fd.write("{:+10.9e}    {:+10.9e}    {:+10.9e}    {:+10.9e}    {:+10.9e}"
        #         .format(time, density[self._item_id], pressure[self._item_id], energy[self._item_id], pseudo[self._item_id]))
        # self._fd.write(os.linesep)
        self._data.append([time, density[self._item_id], pressure[self._item_id], energy[self._item_id], pseudo[self._item_id]])


    @property
    def history_array(self):
        return np.array(self._data)

    def write_fields_history(self):
        """
        write the history data in appropriate file
        """
        # import ipdb ; ipdb.set_trace()
        np.savetxt(self._path, self.history_array, fmt=['%+10.9e', '%+10.9e', '%+10.9e', '%+10.9e', '%+10.9e'],
                   header=self.header)


class NodeTimeData(TimeData):
    """
    Classe définissant un jeu de donnée time history pour un élément
    """

    def __init__(self, id_number, o_path):
        """
        initialisation comme particularisation de TimeData
        :param id_number: numéro du noeud
        :param o_path: chemin pour enregistrer dans le bon dossier
        """
        super(NodeTimeData, self).__init__("node", id_number, o_path)
        self._data = []


    @property
    def header(self):
        """
        crée le fichier pour cell time history
        :return:
        """
        msg = "---Time History for node :" + str(self._item_id) + os.linesep
        msg += "{:^15s}    {:^15s}    {:^15s}".format("Time[s]", "Position[m]", "Vitesse[m/s]")
        msg += os.linesep
        return msg


    def add_time_step_fields(self, time, position, velocity):
        """
        Append the history file with the current cell data
        :var
        - path : str, path to the history data file
        - time : float, real simulation time
        - position : array, position of the node
        - velocity : array, node fields to be added
        - node_number : int, reference of the node of interest to select data
        """
        # self._fd.write("{:+10.9e}    {:+10.9e}    {:+10.9e}"
        #                .format(time, position[self._item_id][0], velocity[self._item_id][0]))
        # self._fd.write(os.linesep)

        self._data.append([time, position[self._item_id], velocity[self._item_id]])

    @property
    def history_array(self):
        return np.array(self._data)

    def write_fields_history(self):
        """
        write the history data in appropriate file
        """
        np.savetxt(self._path, self.history_array, fmt=['%+10.9e', '%+10.9e', '%+10.9e'], header=self.header)

