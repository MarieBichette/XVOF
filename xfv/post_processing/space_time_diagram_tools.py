#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
A class for utilities for plotting a space-time diagram after exploitation of the output database
"""

import pathlib
import numpy as np
import matplotlib.pyplot as plt
from xfv.src.data.data_container import DataContainer
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit


class SpaceTimeDiagramUtilities:
    """
    A class utilities for plotting a space-time diagram
    """

    def __init__(self, path_to_db, verbose):
        """
        Creation of the class
        :param path_to_db : complete path to database hdf5 band (should include all_fields.hdf5)
        :param verbose : to allow a detailed print
        :type str
        :type bool
        """
        self._my_hd = OutputDatabaseExploit(path_to_db)
        self._case = os.path.basename(os.path.dirname(path_to_db))
        self._data_container = DataContainer(os.path.join(os.path.dirname(path_to_db), "XDATA.xml"))
        self._verbose = verbose

    def _data_has_interface(self):
        """
        Return true if an interface exists between node and projectile
        :return: bool
        """
        try:
            interface_position = self.data_container.geometric.initial_interface_position
            return True
        except ValueError:
            return False

    def compute_interface_cell_id(self):
        """
        Return the cell id and the position at the interface between projectile and target
        (last cell in the projectile)
        :return: tuple(cell_interface_id (int), interface_position(float)
        """
        interface_position = self.data_container.geometric.initial_interface_position
        if self._verbose:
            print("Interface position : " + str(interface_position))

        # Look the index of cell of the projectile / target interface
        t = 0
        node_coord = self.my_hd.extract_field_at_time("NodeCoordinates", t)[:-1].flatten()
        cell_size = self.my_hd.extract_field_at_time("CellSize", t)[:]
        cell_coord = (node_coord + cell_size / 2.)
        indice_interface = 0
        for i_cell in range(cell_size.shape[0]):
            if cell_coord[i_cell] <= interface_position:
                indice_interface += 1
        return indice_interface, interface_position

    def plot_interface(self, coordinates_array, time_array, i_fig=1):
        """
        Plot the interface between interface and projectile
        :param coordinates_array: coordinates
        :param time_array: time
        :param i_fig : figure index to add the interface line to
        :return:
        """
        index, position = self.compute_interface_cell_id()
        if self._verbose:
            print("Put the interface between projectile and target at " + str(position) +
                  "(= cell " + str(indice) + ")")
        plt.figure(i_fig)
        plt.plot(coordinates_array[:, index], time_array[:, index], color='black')

    def build_xyz_map_for_contourf_plot(self, field_type):
        """
        Build a map to plot the space time diagram with contourf matplotlib function
        :param field_type : string for field type of interest.
        Should match a key of dictionary OutputDatabaseExploit.field_type_converter
        :type field_type : str
        :return: tuple(X, Y, Z) map
        X = coordinates, Y = time et Z = field to map
        """
        classical_fname, additional_fname = OutputDatabaseExploit.field_type_converter[field_type]

        # 1) Initialiser les maps
        # astuce: on crée plusieurs map pour chaque "partie" de la barre rompue et on trace
        # pour les cells pas encore rompues, on dédouble la valeur
        # -> on se retrouve après construction des tableaux avec
        # dim_x qui vaut nbr_cell + nbr_disc à la fin du calcul

        # 1.a) Récupérer les disc finales
        final_cell_status = \
            self.my_hd.extract_field_at_time("CellStatus", self.my_hd.saved_times[-1])
        final_ruptured_cell_id = np.where(final_cell_status)[0]
        final_ruptured_cell_id = np.sort(final_ruptured_cell_id)

        # 1.b) calculer les tailles des tableaux
        dim_x_0 = self.my_hd.extract_field_at_time(classical_fname, 0.).shape[0]
        dim_x_final = dim_x_0 + len(final_ruptured_cell_id)
        dim_y = self.my_hd.nb_saved_times

        # 1.c) Créer les tableaux
        X = np.zeros([dim_y, dim_x_final])
        Y = np.zeros([dim_y, dim_x_final])
        Z = np.zeros([dim_y, dim_x_final])

        # 2) Build line after line the color map X, Y, Z
        for i_temps in range(self.my_hd.nb_saved_times):
            t = self.my_hd.saved_times[i_temps]
            x, y, z = self._build_classical_xyz_map_for_contourf_plot_at_time(t, dim_x_0)
            self._include_enrichment(x, y, z, t, final_ruptured_cell_id)  # modifies x, y, z
            X[i_temps, :] = x * 1.e3
            Y[i_temps, :] = y * 1.e6
            Z[i_temps, :] = z
        return X, Y, Z

    def plot_section_of_color_map(self, X, Y, Z, begin=None, end=None,
                                  n_colors= 500, vmin=0, vmax=1):
        """
        Plot a line of the contourf map with defined boundaries in space begin and end
        :param X: coordinates array
        :param Y: time array
        :param Z: field array
        :param begin: left boundary id
        :param end: right boundary id
        :param n_colors: param of the plot : number of color in the colorbar
        :param vmin: param of the plot : minimum of the colorbar
        :param vmax: param of the plot : maximum of the colorbar
        :return:
        """
        if begin is not None and end is not None:
            plt.contourf(X[:, begin:end + 1], Y[:, begin:end + 1], Z[:, begin:end + 1],
                         n_colors, vmin=vmin, vmax=vmax)
            if self._verbose:
                print("Plot from " + str(begin) + " to " + str(end) + "inclus")

        elif begin is None and end is not None:
            plt.contourf(X[:, :end + 1], Y[:, :end + 1], Z[:, :end + 1],
                         n_colors, vmin=vmin, vmax=vmax)
            if self._verbose:
                print("Plot from the beginning to " + str(end) + "inclus")

        elif begin is not None and end is None:
            plt.contourf(X[:, begin:], Y[:, begin:], Z[:, begin:], n_colors, vmin=vmin, vmax=vmax)
            if self._verbose:
                print("Plot from " + str(begin) + " to the end")

        elif begin is None and end is None:
            plt.contourf(X[:, :], Y[:, :], Z[:, :], n_colors, vmin=vmin, vmax=vmax)
            if self._verbose:
                print("Plot all data")

    def _build_classical_xyz_map_for_contourf_plot_at_time(self, t, dim_x):
        """
        Extract information out of the database at time t
        and return a convenient storage for contour_f
        :param t: time
        :param dim_x :
        :type: float
        :type : int
        :return: tuple(array(coordinates), array_ones * time, array(field))
        """
        node_coord = self.my_hd.extract_field_at_time("NodeCoordinates", t)[:-1].flatten()
        cell_size = self.my_hd.extract_field_at_time("CellSize", t)[:]
        x = (node_coord + cell_size / 2.)
        y = np.ones([dim_x]) * t
        z = self.my_hd.extract_field_at_time(classical_fname, t)[:]
        return x, y, z

    def _include_enrichment(self, x, y, z, t, final_ruptured_cell_id):
        """
        Modifies the x,y,z array to include enrichment either by dedoubling the data if the disc is
        not created at time t yet, or by inserting the disc data if the disc exists at time t
        :param x:
        :param y:
        :param z:
        :param t:
        :return:
        """
        # Initialisation of useful data
        offset = 0
        count_active_disc = 0
        current_cell_status = self.my_hd.extract_field_at_time("CellStatus", t)[:]

        if current_cell_status.any():
            # Mutualisation de variables qui seront demandées item par item plus tard
            left_size = self.my_hd.extract_field_at_time("AdditionalLeftSize", t)[:]
            right_size = self.my_hd.extract_field_at_time("AdditionalRightSize", t)[:]
            right_field = self.my_hd.extract_field_at_time(additional_fname, t)[:]
            left_size = np.sort(left_size, 0)
            right_size = np.sort(right_size, 0)
            right_field = np.sort(right_field, 0)

        # Boucle on all disc that will be created at final time
        for i_cell_index in final_ruptured_cell_id:
            moving_index = int(i_cell_index) + offset
            if not current_cell_status[i_cell_index]:
                # la cell n'est pas encore enrichie mais va l'être dans la suite  du calcul :
                # on dédouble la valeur classique en prévision de l'insertion des
                x = np.insert(x, moving_index + 1, x[moving_index])
                y = np.insert(y, moving_index + 1, y[moving_index])
                z = np.insert(z, moving_index + 1, z[moving_index])

            else:
                # si la cell est enrichie : on insère les données pour la partie droite :
                # * pour x : on modifie la taille de la cell de gauche + on insère la cell de droite
                left_size_for_cell_i = left_size[count_active_disc, 1]
                right_size_for_cell_i = right_size[count_active_disc, 1]
                x[moving_index] += left_size_for_cell_i / 2. - cell_size[i_cell_index] / 2.
                right_coordinate = x[moving_index + 1] - right_size_for_cell_i / 2. - cell_size[
                    i_cell_index + 1] / 2.
                x = np.insert(x, moving_index + 1, [right_coordinate])

                # * pour y : on insère une nouvelle valeur de t
                y = np.insert(y, moving_index + 1, [t])

                # * pour z : on insère le champ droite
                right_field_for_cell_i = right_field[count_active_disc, 1]
                z = np.insert(z, moving_index + 1, [right_field_for_cell_i])
                count_active_disc += 1
            offset += 1
