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
        self._data_container = DataContainer(os.path.join(os.path.dirname(path_to_db),
                                                          "XDATA.json"))
        self._verbose = verbose
        self.first_enr_time = dict()

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

    def build_xyz_map_for_contourf_plot(self, field_type: str, final_ruptured_cell_id):
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

        # Create arrays as lists
        coord_array = []
        time_array = []
        field_array = []

        # 2) Build line after line the color map X, Y, Z
        for i_temps in range(self._my_hd.nb_saved_times):
            time = self._my_hd.saved_times[i_temps]
            cell_size = self._my_hd.extract_field_at_time("CellSize", time)[:]
            coord_array_at_time, time_array_at_time, field_array_at_time = \
                self._build_classical_xyz_map_for_contourf_plot_at_time(time, classical_fname)
            self._include_enrichment(coord_array_at_time, time_array_at_time, field_array_at_time,
                                     time, cell_size, final_ruptured_cell_id, additional_fname)
            scaled_coord = [the_coord * 1.e3 for the_coord in coord_array_at_time]
            scaled_time = [the_time * 1.e6 for the_time in time_array_at_time]
            coord_array.append(scaled_coord)
            time_array.append(scaled_time)
            field_array.append(field_array_at_time)
        return coord_array, time_array, field_array  # X, Y, Z

    def plot_section_of_color_map(self, coord, time, field, options):
        """
        Plot a line of the contourf map with defined boundaries in space begin and end
        :param coord: coordinates array
        :param time: time array
        :param field: field array
        :param options : tuple(int, float, float) with :
                - n_colors: param of the plot : number of color in the color bar
                - field_min: param of the plot : minimum of the color bar
                - field_max: param of the plot : maximum of the color bar
        :return:
        """
        plt.contourf(coord, time, field, options.n_colors,
                     vmin=options.field_min, vmax=options.field_max)

    def plot_geometry_boundaries(self, coord, time):
        """
        Plot the discontinuities boundary in black line
        :param coord :  coord array of the contourf function
        :param time : time array of the contourf function
        """
        plt.plot(coord[:, 0], time[:, 0], color='black', linewidth=0.5)  # contour bord g
        plt.plot(coord[:, -1], time[:, -1], color='black', linewidth=0.5)  # contour bord d

    def plot_disc_boundaries(self, coord, time, cell_id_before_offset, cell_id_after_offset):
        """
        Plot the discontinuities boundary in black line
        :param coord :  coord array of the contourf function
        :param time : time array of the contourf function
        :param cell_id_before_offset : cell id of the current disc to plot
        :param cell_id_after_offset : take into account the offset of previous discs
        """
        disc_time = self.first_enr_time[cell_id_before_offset]
        index_time = self.get_enrichment_time_index(disc_time)
        n_id = cell_id_after_offset + 1  # next id
        plt.plot(coord[index_time:, cell_id_after_offset],
                 time[index_time:, cell_id_after_offset],
                 color='black', linewidth=0.5)  # contour bord g
        plt.plot(coord[index_time:, n_id], time[index_time:, n_id],
                 color='black', linewidth=0.5)  # contour bord d

    def _build_classical_xyz_map_for_contourf_plot_at_time(self, time: float,
                                                           classical_fname: str):
        """
        Extract information out of the database at time t
        and return a convenient storage for contour_f
        :param time: time
        :param classical_fname : name of the classical field (starting by Classical...)
        :return: tuple(array(coordinates), array_ones * time, array(field))
        """
        node_coord = self._my_hd.extract_field_at_time("NodeCoordinates", time)[:-1].flatten()
        cell_size = self._my_hd.extract_field_at_time("CellSize", time)[:]
        coord_array = (node_coord + cell_size / 2.).tolist()
        time_array = (np.ones_like(cell_size) * time).tolist()
        field_array = (self._my_hd.extract_field_at_time(classical_fname, time)[:]).tolist()
        return coord_array, time_array, field_array

    def _include_enrichment(self, coord_array: list, time_array: list, field_array: list, time,
                            cell_size, final_ruptured_cell_id: np.array, additional_fname: str):
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
        offset: int = 0
        count_active_disc: int = 0
        current_cell_status = self._my_hd.extract_field_at_time("CellStatus", time)[:]

        if current_cell_status.any():
            # Initialisation of variables which will be used later in item loop
            left_size = self._my_hd.extract_field_at_time("AdditionalLeftSize", time)[:]
            right_size = self._my_hd.extract_field_at_time("AdditionalRightSize", time)[:]
            right_field = self._my_hd.extract_field_at_time(additional_fname, time)[:]
            sorted_disc = np.argsort(left_size[:, 0])
            left_size = left_size[sorted_disc, 1].tolist()
            right_size = right_size[sorted_disc, 1].tolist()
            right_field = right_field[sorted_disc, 1].tolist()
        else:
            left_size, right_size = None, None
            right_field = None

        # Loop on all disc that will be created at final time
        for i_cell_index in final_ruptured_cell_id:  # ids of cells that are cracked at the end
            moving_index = int(i_cell_index) + offset
            if not current_cell_status[i_cell_index]:
                # Current_cell is not enriched yet but will be in the future
                # Double the value of classical field to anticipate the future fracture
                coord_array.insert(moving_index + 1, coord_array[moving_index])
                time_array.insert(moving_index + 1, time_array[moving_index])
                field_array.insert(moving_index + 1, field_array[moving_index])

            else:
                if i_cell_index not in self.first_enr_time:
                    self.first_enr_time[i_cell_index] = time
                # If the cell is enriched, insert the value for right part of the cracked cell
                # * for coord_array : add left_size and insert right size of cell
                left_size_for_cell_i = left_size[count_active_disc]
                right_size_for_cell_i = right_size[count_active_disc]
                # Here, we compute the coordinates of the left boundary of the discontinuity
                # instead of the center of the left part of cracked cell for representativeness of
                # the diagram
                coord_array[moving_index] += \
                    left_size_for_cell_i - cell_size[i_cell_index] / 2.
                # In the same idea, we compute the coordinates of the right boundary of the
                # discontinuity instead of the center of the right part of cracked cell for
                # representativeness of the diagram
                right_coordinate = (coord_array[moving_index + 1]
                                    - right_size_for_cell_i - cell_size[i_cell_index + 1] / 2.)
                coord_array.insert(moving_index + 1, right_coordinate)

                # * for time_array : insert timed_right_
                time_array.insert(moving_index + 1, time)

                # * for field_array : insert value of right field
                right_field_for_cell_i = right_field[count_active_disc]
                field_array.insert(moving_index + 1, right_field_for_cell_i)
                count_active_disc += 1
            offset += 1

    def get_enrichment_time_index(self, time):
        """
        Returns the index of first time with discontinuities
        :param time : time to search
        """
        return self._my_hd.saved_times.index(time)
