#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
A class for utilities for plotting a space-time diagram after exploitation of the output database
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from xfv.src.data.data_container import DataContainer
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit


class SpaceTimeDiagramTools:
    """
    A class utilities for plotting a space-time diagram
    """

    def __init__(self, path_to_db: str, verbose: bool):
        """
        Creation of the class
        :param path_to_db : complete path to database hdf5 band (should include all_fields.hdf5)
        :param verbose : to allow a detailed print
        """
        self._my_hd = OutputDatabaseExploit(path_to_db)
        self._case = os.path.basename(os.path.dirname(path_to_db))
        self._data_container = DataContainer(os.path.join(os.path.dirname(path_to_db), "XDATA.json"))
        self._verbose = verbose

    def data_has_interface(self):
        """
        Return true if an interface exists between node and projectile
        :return: bool
        """
        try:
            interface_position = self._data_container.geometric.initial_interface_position
            return True
        except ValueError:
            return False

    def compute_interface_cell_id(self):
        """
        Return the cell id and the position at the interface between projectile and target
        (last cell in the projectile)
        :return: tuple(cell_interface_id (int), interface_position(float)
        """
        interface_position = self._data_container.geometric.initial_interface_position
        if self._verbose:
            print("Interface position : " + str(interface_position))

        # Look the index of cell of the projectile / target interface
        time = 0
        node_coord = self._my_hd.extract_field_at_time("NodeCoordinates", time)[:-1].flatten()
        cell_size = self._my_hd.extract_field_at_time("CellSize", time)[:]
        cell_coord = (node_coord + cell_size / 2.)
        indice_interface = 0
        for i_cell in range(cell_size.shape[0]):
            if cell_coord[i_cell] <= interface_position:
                indice_interface += 1
        return indice_interface, interface_position

    def plot_interface(self, coordinates_array, time_array):
        """
        Plot the interface between interface and projectile
        :param coordinates_array: coordinates
        :param time_array: time
        :return:
        """
        # In theory, interface is a node but the array of coordinates is based on cell
        index, position = self.compute_interface_cell_id()
        if self._verbose:
            print("Put the interface between projectile and target at " + str(position) +
                  "(= cell " + str(index) + ")")
        plt.plot(coordinates_array[:, index], time_array[:, index], color='black')

    def build_xyz_map_for_contourf_plot(self, field_type: str, final_ruptured_cell_id: int):
        """
        Build a map to plot the space time diagram with contourf matplotlib function
        :param field_type : string for field type of interest.
        Should match a key of dictionary OutputDatabaseExploit.field_type_converter
        :param final_ruptured_cell_id : array of the enriched cells at final time
        :return: tuple(X, Y, Z) map
        X = coordinates, Y = time et Z = field to map
        """
        classical_fname, additional_fname = OutputDatabaseExploit.field_type_converter[field_type]

        # 1) Initialize maps
        # trick: create several data for each part of the cracked bar
        # For the cells to be cracked (but not yet cracked), double value of classical field
        # At the end, build array with dimension dim_x = nbr_cell + nbr_disc at the end of
        # calculations
        # 1.a) Compute array size
        dim_x_0 = self._my_hd.extract_field_at_time(classical_fname, 0.).shape[0]
        dim_x_final = dim_x_0 + len(final_ruptured_cell_id)
        dim_y = self._my_hd.nb_saved_times

        # 1.b) Create arrays
        coord_array = np.zeros([dim_y, dim_x_final])
        time_array = np.zeros([dim_y, dim_x_final])
        field_array = np.zeros([dim_y, dim_x_final])

        # 2) Build line after line the color map X, Y, Z
        for i_temps in range(self._my_hd.nb_saved_times):
            time = self._my_hd.saved_times[i_temps]
            cell_size = self._my_hd.extract_field_at_time("CellSize", time)[:]
            coord_array_at_time, time_array_at_time, field_array_at_time = \
                self._build_classical_xyz_map_for_contourf_plot_at_time(time, dim_x_0,
                                                                        classical_fname)
            modified_coord_array_t, modified_time_array_t, modified_field_array_t = \
                self._include_enrichment(coord_array_at_time, time_array_at_time,
                                         field_array_at_time, time, cell_size,
                                         final_ruptured_cell_id,
                                         additional_fname)  # modifies coord_at_time, y, z
            coord_array[i_temps, :] = modified_coord_array_t * 1.e3
            time_array[i_temps, :] = modified_time_array_t * 1.e6
            field_array[i_temps, :] = modified_field_array_t
        return coord_array, time_array, field_array  # X, Y, Z

    def plot_section_of_color_map(self, coord, time, field, options, begin=None, end=None):
        """
        Plot a line of the contourf map with defined boundaries in space begin and end
        :param coord: coordinates array
        :param time: time array
        :param field: field array
        :param begin: left boundary id
        :param end: right boundary id
        :param options : tuple(int, float, float) with :
                - n_colors: param of the plot : number of color in the color bar
                - field_min: param of the plot : minimum of the color bar
                - field_max: param of the plot : maximum of the color bar
        :return:
        """
        if begin is not None and end is not None:
            plt.contourf(coord[:, begin:end + 1], time[:, begin:end + 1], field[:, begin:end + 1],
                         options.n_colors, vmin=options.field_min, vmax=options.field_max)
            plt.plot(coord[:, begin], time[:, begin], color='black', linewidth=0.5)  # contour g
            plt.plot(coord[:, end], time[:, end], color='black', linewidth=0.5)  # contour d
            if self._verbose:
                print("Plot from " + str(begin) + " to " + str(end) + "inclus")

        elif begin is None and end is not None:
            plt.contourf(coord[:, :end + 1], time[:, :end + 1], field[:, :end + 1],
                         options.n_colors, vmin=options.field_min, vmax=options.field_max)
            plt.plot(coord[:, 0], time[:, 0], color='black', linewidth=0.5)  # contour g
            plt.plot(coord[:, end], time[:, end], color='black', linewidth=0.5)  # contour d
            if self._verbose:
                print("Plot from the beginning to " + str(end) + "inclus")

        elif begin is not None and end is None:
            plt.contourf(coord[:, begin:], time[:, begin:], field[:, begin:],
                         options.n_colors, vmin=options.field_min, vmax=options.field_max)
            plt.plot(coord[:, begin], time[:, begin], color='black', linewidth=0.5)  # contour g
            plt.plot(coord[:, -1], time[:, -1], color='black', linewidth=0.5)  # contour d
            if self._verbose:
                print("Plot from " + str(begin) + " to the end")

        elif begin is None and end is None:
            plt.contourf(coord[:, :], time[:, :], field[:, :],
                         options.n_colors, vmin=options.field_min, vmax=options.field_max)
            plt.plot(coord[:, 0], time[:, 0], color='black', linewidth=0.5)  # contour g
            plt.plot(coord[:, -1], time[:, -1], color='black', linewidth=0.5)  # contour d
            # no contour
            if self._verbose:
                print("Plot data for all geometry")

    def _build_classical_xyz_map_for_contourf_plot_at_time(self, time: float, dim_x: int,
                                                           classical_fname: str):
        """
        Extract information out of the database at time t
        and return a convenient storage for contour_f
        :param time: time
        :param dim_x : nbr of cells
        :param classical_fname : name of the classical field (starting by Classical...)
        :return: tuple(array(coordinates), array_ones * time, array(field))
        """
        node_coord = self._my_hd.extract_field_at_time("NodeCoordinates", time)[:-1].flatten()
        cell_size = self._my_hd.extract_field_at_time("CellSize", time)[:]
        coord_array = (node_coord + cell_size / 2.)
        time_array = np.ones([dim_x]) * time
        field_array = self._my_hd.extract_field_at_time(classical_fname, time)[:]
        return coord_array, time_array, field_array

    def _include_enrichment(self, coord_array, time_array, field_array, time, cell_size,
                            final_ruptured_cell_id: np.array, additional_fname: str):
        """
        Modifies the x,y,z array to include enrichment either by dedoubling the data if the disc is
        not created at time t yet, or by inserting the disc data if the disc exists at time t
        :param coord_array: coord array for current time
        :param time_array: time array for current time
        :param field_array: field array for current time
        :param time: time
        :param cell_size : field of cell classical size
        :param final_ruptured_cell_id : id of the enriched cell after offset
        :param additional_fname : name of the enriched field (starting with Additional...)
        :return:
        """
        # Initialisation of useful data
        offset = 0
        count_active_disc = 0
        current_cell_status = self._my_hd.extract_field_at_time("CellStatus", time)[:]

        modified_coord_array = np.copy(coord_array)
        modified_time_array = np.copy(time_array)
        modified_field_array = np.copy(field_array)

        if current_cell_status.any():
            # Initialisation of variables which will be used later in item loop
            left_size = self._my_hd.extract_field_at_time("AdditionalLeftSize", time)[:]
            right_size = self._my_hd.extract_field_at_time("AdditionalRightSize", time)[:]
            right_field = self._my_hd.extract_field_at_time(additional_fname, time)[:]
            left_size = np.sort(left_size, 0)
            right_size = np.sort(right_size, 0)
            right_field = np.sort(right_field, 0)
        else:
            left_size, right_size = None, None
            right_field = None, None

        # Loop on all disc that will be created at final time
        for i_cell_index in final_ruptured_cell_id:
            moving_index = int(i_cell_index) + offset
            if not current_cell_status[i_cell_index]:
                # Current_cell is not enriched yet but will be in the future
                # Double the value of classical field to anticipate the future fracture
                modified_coord_array = np.insert(modified_coord_array, moving_index + 1,
                                                 modified_coord_array[moving_index])
                modified_time_array = np.insert(modified_time_array, moving_index + 1,
                                                modified_time_array[moving_index])
                modified_field_array = np.insert(modified_field_array, moving_index + 1,
                                                 modified_field_array[moving_index])

            else:
                # If the cell is enriched, insert the value for right part of the cracked cell
                # * for coord_array : add left_size and insert right size of cell
                left_size_for_cell_i = left_size[count_active_disc, 1]
                right_size_for_cell_i = right_size[count_active_disc, 1]
                # Here, we compute the coordinates of the left boundary of the discontinuity
                # instead of the center of the left part of cracked cell for representativeness of
                # the diagram
                modified_coord_array[moving_index] += \
                    left_size_for_cell_i - cell_size[i_cell_index] / 2.
                # In the same idea, we compute the coordinates of the right boundary of the
                # discontinuity instead of the center of the right part of cracked cell for
                # representativeness of the diagram
                right_coordinate = modified_coord_array[moving_index + 1] \
                                   - right_size_for_cell_i - cell_size[i_cell_index + 1] / 2.
                modified_coord_array = np.insert(modified_coord_array, moving_index + 1, [right_coordinate])

                # * for time_array : insert time
                modified_time_array = np.insert(modified_time_array, moving_index + 1, [time])

                # * for field_array : insert value of right field
                right_field_for_cell_i = right_field[count_active_disc, 1]
                modified_field_array = np.insert(modified_field_array, moving_index + 1,
                                                 [right_field_for_cell_i])
                count_active_disc += 1
            offset += 1

        return modified_coord_array, modified_time_array, modified_field_array
