#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
A class for utilities for plotting a march diagram after exploitation of the output database
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from xfv.src.data.data_container import DataContainer
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit


class SpaceTimeDiagramUtilities:

    def __init__(self, case):
        """

        """
        hdf5_file = os.path.join(case.directory_name, "all_fields.hdf5")
        self.my_hd = OutputDatabaseExploit(hdf5_file)
        self.case = case
        self.data_container = DataContainer(os.path.join(self.case.directory_name, "XDATA.xml"))

    def data_has_interface(self):
        """
        Return the cell id at the interface between projectile and target (last cell in the projectile)
        :return: cell_interface_id (int)
        """
        try:
            interface_position = self.data_container.geometric.initial_interface_position
            return True
        except ValueError:
            return False

    def compute_interface_cell_id(self):
        """
        Return the cell id at the interface between projectile and target (last cell in the projectile)
        :return: cell_interface_id (int)
        """
        interface_position = self.data_container.geometric.initial_interface_position

        # Analyse valable à t=0 (pas encore d'enrichissement)
        t = 0
        node_coord = self.my_hd.extract_field_at_time("NodeCoordinates", t)[:-1].flatten()
        cell_size = self.my_hd.extract_field_at_time("CellSize", t)[:]
        x = (node_coord + cell_size / 2.)

        # On cherche l'indice qui correspond à l'interface
        indice_interface = 0
        for i_cell in range(cell_size.shape[0]):
            if x[i_cell] <= interface_position:
                indice_interface += 1

        return indice_interface, interface_position

    def plot_interface(self, X, Y, i_fig=1):
        """
        Trace sur la figure i_fig la position de l'interface projectile / cible
        :param X:
        :param Y:
        :param i_fig : indice de la figure
        :return:
        """
        indice, position = self.compute_interface_cell_id()
        print("Positionnement de l'interface projectile / cible à " + str(position) + "(= cell " + str(indice) + ")")
        plt.figure(i_fig)
        plt.plot(X[:, indice], Y[:, indice], color='black')


    def build_XYZ_map_for_contourf_plot(self, field_type):
        """
        :param field_type : string forfield type of interest.
        Should match a keyof dictionary OutputDatabaseExploit.field_type_converter
        :return: X Y Z map X = positions, Y = temps et Z = champ à tracer
        """

        classical_fname, additional_fname = OutputDatabaseExploit.field_type_converter[field_type]

        final_cell_status = self.my_hd.extract_field_at_time("CellStatus", self.my_hd.saved_times[-1])
        final_ruptured_cell_id = np.where(final_cell_status)[0]
        final_ruptured_cell_id = np.sort(final_ruptured_cell_id)

        dim_x = self.my_hd.extract_field_at_time(classical_fname, 0.).shape[0]
        dim_x_final = dim_x + len(final_ruptured_cell_id)
        dim_y = self.my_hd.nb_saved_times

        X = np.zeros([dim_y, dim_x_final])
        Y = np.zeros([dim_y, dim_x_final])
        Z = np.zeros([dim_y, dim_x_final])
        # astuce: on crée plusieurs map pour chaque "partie" de la barre rompue et on trace
        # pour les cells pas encore rompues, on dédouble la valeur
        # -> on se retrouve après construction des tableaux avec dim_x qui vaut nbr_cell + nbr_disc à la fin du calcul

        for i_temps in range(self.my_hd.nb_saved_times):
            t = self.my_hd.saved_times[i_temps]
            offset = 0
            compteur_disc_active = 0
            node_coord = self.my_hd.extract_field_at_time("NodeCoordinates", t)[:-1].flatten()
            cell_size = self.my_hd.extract_field_at_time("CellSize", t)[:]
            x = (node_coord + cell_size / 2.)
            y = np.ones([dim_x]) * t
            z = self.my_hd.extract_field_at_time(classical_fname, t)[:]

            current_cell_status = self.my_hd.extract_field_at_time("CellStatus", t)[:]

            if current_cell_status.any():
                left_size = self.my_hd.extract_field_at_time("AdditionalLeftSize", t)[:]
                right_size = self.my_hd.extract_field_at_time("AdditionalRightSize", t)[:]
                right_field = self.my_hd.extract_field_at_time(additional_fname, t)[:]
                left_size = np.sort(left_size, 0)
                right_size = np.sort(right_size, 0)
                right_field = np.sort(right_field, 0)

            for i_cell_index in final_ruptured_cell_id:
                moving_index = int(i_cell_index) + offset
                if not current_cell_status[i_cell_index]:
                    # la cell n'est pas encore enrichie mais va l'être dans la suite  du calcul :
                    # on dédouble la valeur classique en prévision de l'insertion des
                    x = np.insert(x, moving_index + 1, x[moving_index])
                    y = np.insert(y, moving_index + 1, y[moving_index])
                    z = np.insert(z, moving_index + 1, z[moving_index])

                else:
                    # si la cell est enrichie : on insère la partie droite :
                    # pour x : on modifie la taille de la cell de gauche et on insère la cell de droite
                    left_size_for_cell_i = left_size[compteur_disc_active, 1]
                    right_size_for_cell_i = right_size[compteur_disc_active, 1]
                    x[moving_index] += left_size_for_cell_i / 2. - cell_size[i_cell_index] / 2.
                    right_coordinate = x[moving_index + 1] - right_size_for_cell_i / 2. - cell_size[i_cell_index + 1] / 2.
                    x = np.insert(x, moving_index + 1, [right_coordinate])

                    # pour y : on insère une nouvelle valeur de t (on s'en fiche puisque y est un vecteur de constantes)
                    y = np.insert(y, moving_index + 1, [t])
                    # pour z : on insère le champ droite
                    right_field_for_cell_i = right_field[compteur_disc_active, 1]
                    z = np.insert(z, moving_index + 1, [right_field_for_cell_i])
                    compteur_disc_active += 1
                offset += 1

            X[i_temps, :] = x * 1.e3
            Y[i_temps, :] = y * 1.e6
            Z[i_temps, :] = z
        return X, Y, Z

    def plot_section_of_color_map(self, X, Y, Z, begin=None, end=None, fig=1, n_colors= 500, vmin=0, vmax=1):

        plt.figure(fig)

        if begin is not None and end is not None:
            plt.contourf(X[:, begin:end + 1], Y[:, begin:end + 1], Z[:, begin:end + 1], n_colors, vmin=vmin, vmax=vmax)
            print("Plot from " + str(begin) + " to " + str(end) + "inclus")

        elif begin is None and end is not None:
            plt.contourf(X[:, :end + 1], Y[:, :end + 1], Z[:, :end + 1], n_colors, vmin=vmin, vmax=vmax)
            print("Plot from the beginning to " + str(end) + "inclus")

        elif begin is not None and end is None:
            plt.contourf(X[:, begin:], Y[:, begin:], Z[:, begin:], n_colors, vmin=vmin, vmax=vmax)
            print("Plot from " + str(begin) + " to the end")

        elif begin is None and end is None:
            plt.contourf(X[:, :], Y[:, :], Z[:, :], n_colors, vmin=vmin, vmax=vmax)
            print("Plot all data")



