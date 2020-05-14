# -*- coding: utf-8 -*-
"""
A script plotting a space-time diagram after exploitation of the output database
"""

import argparse
import pathlib
from collections import namedtuple
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit
from xfv.post_processing.tools.space_time_diagram_tools import SpaceTimeDiagramTools


def create_figure():
    """
    Creation of the color map figure
    """
    fig = plt.figure(1)
    plt.xlabel("Coordinates [mm]", fontsize=18)
    plt.ylabel("Time [mu s]", fontsize=18)
    if ARGS.gradient:
        the_title = "Space time {:} gradient diagram".format(ARGS.field)
        # Definition of a color map for gradients plot
        color_dict = {'red': ((0., 0., 0.), (0.5, 0.75, 0.75), (1., 1., 1.)),
                      'green': ((0., 0., 0.), (0.5, 0.75, 0.75), (1., 0., 0.)),
                      'blue': ((0., 0.5, 0.5), (0.5, 0.75, 0.75), (1., 0., 0.))}
        my_cmap = LinearSegmentedColormap('custom_cmap', color_dict)
        plt.register_cmap(cmap=my_cmap)
        plt.set_cmap('custom_cmap')  # set the colormap just created
    else:
        the_title = "Diagram {:}".format(ARGS.case)
    plt.title(the_title, fontweight='bold')
    return fig


def show_figure(fig, min_field: float, max_field: float):
    """
    Show figure
    :param fig : figure
    :param min_field: lower bound of the color bar
    :param max_field: upper bound of the color bar
    """
    # Color bar legend
    axes, _ = matplotlib.colorbar.make_axes(fig.gca())
    color_bar = matplotlib.colorbar.ColorbarBase(
        axes, norm=matplotlib.colors.Normalize(vmin=min_field, vmax=max_field))
    if ARGS.gradient:
        color_bar.set_label("{:} gradient".format(ARGS.field), fontsize=18)
    else:
        color_bar.set_label(ARGS.field, fontsize=18)
    # Show plot
    plt.legend()
    plt.show()


def run():
    """
    Run post processing
    """
    path_to_db = pathlib.Path.cwd().joinpath("..", "tests", ARGS.case, ARGS.output_filename)
    my_hd = OutputDatabaseExploit(path_to_db)

    # ----------------------------------------------------------------
    # Get the final number of created discontinuities
    # ----------------------------------------------------------------
    final_cell_status = my_hd.extract_field_at_time("CellStatus", my_hd.saved_times[-1])
    ruptured_cell_id_before_offset = np.where(final_cell_status)[0]
    ruptured_cell_id_before_offset = np.sort(ruptured_cell_id_before_offset)
    # ----------------------------------------------------------------
    # Get data to plot the diagram
    # ----------------------------------------------------------------
    diagram_tools = SpaceTimeDiagramTools(path_to_db, ARGS.verbose)
    coord_list, time_list, field_list = \
        diagram_tools.build_xyz_map_for_contourf_plot(ARGS.field, ruptured_cell_id_before_offset)
    coord_array = np.array(coord_list)
    time_array = np.array(time_list)
    field_array = np.array(field_list)
    min_field = np.min(field_array)
    max_field = np.max(field_array)

    if ARGS.gradient:
        field_array = np.gradient(field_array, axis=0)  # time derivative of the field (axis = 0)
        min_field = -5.e8
        max_field = 5.e8
        if ARGS.verbose:
            print("Scaling the color map between {:} and {:}".format(min_field, max_field))

    PlotOptions = namedtuple("PlotOptions", ["n_colors", "field_min", "field_max"])
    plot_options = PlotOptions(N_COLORS, min_field, max_field)

    # ----------------------------------------------------------------
    # Tracé color map 2D
    # ----------------------------------------------------------------
    if ARGS.verbose:
        print("Plot the color map")

    # Tracé jusqu'au temps d'appartion des disc :
    first_time = min([diagram_tools.first_enr_time[k] for k in diagram_tools.first_enr_time])
    time_index = diagram_tools.get_enrichment_time_index(first_time)
    print(time_index)
    diagram_tools.plot_section_of_color_map(coord_array[:time_index + 1, :],
                                            time_array[:time_index + 1, :],
                                            field_array[:time_index + 1, :],
                                            plot_options)

    # Apparition des disc :
    if len(ruptured_cell_id_before_offset) >= 1:
        # Offset of the cracked cell ids to be conservative with the number of items of final arrays
        ruptured_cell_id_after_offset = \
            ruptured_cell_id_before_offset + list(range(0, len(ruptured_cell_id_before_offset)))
        if ARGS.verbose:
            print("List of cracked cells :" + str(ruptured_cell_id_before_offset))
            print("=> List of cracked cells after offset:" + str(ruptured_cell_id_after_offset))

        # First part from 0 to first discontinuity
        end = ruptured_cell_id_after_offset[0]  # left boundary of first discontinuity
        diagram_tools.plot_section_of_color_map(coord_array[time_index:, :end + 1],
                                                time_array[time_index:, :end + 1],
                                                field_array[time_index:, :end + 1],
                                                plot_options)
        if ARGS.verbose:
            print("Plot from the beginning to " + str(end) + "included")
        offset = 1

        # From one discontinuity to another (plot from time index where one discontinuity exists)
        for i in range(0, len(ruptured_cell_id_before_offset) - 1):  # boucle sur les disc
            i_rupture_index_before_offset = ruptured_cell_id_before_offset[i]
            i_rupture_index_after_offset = ruptured_cell_id_after_offset[i]
            begin = i_rupture_index_after_offset + 1  # right boundary of the current disc
            end = ruptured_cell_id_after_offset[offset]  # left boundary of the next disc
            diagram_tools.plot_section_of_color_map(coord_array[time_index:, begin:end + 1],
                                                    time_array[time_index:, begin:end + 1],
                                                    field_array[time_index:, begin:end + 1],
                                                    plot_options)
            if ARGS.verbose:
                print("Plot from " + str(begin) + " to " + str(end) + "included")

            diagram_tools.plot_disc_boundaries(coord_array, time_array,
                                               i_rupture_index_before_offset,
                                               i_rupture_index_after_offset)
            offset += 1

        # From the last discontinuity to the end
        begin = ruptured_cell_id_after_offset[-1] + 1  # left boundary of the last disc
        diagram_tools.plot_section_of_color_map(coord_array[time_index:, begin:],
                                                time_array[time_index:, begin:],
                                                field_array[time_index:, begin:],
                                                plot_options)
        if ARGS.verbose:
            print("Plot from " + str(begin) + " to the end")

    else:
        # Simple plot with no discontinuities
        diagram_tools.plot_section_of_color_map(coord_array, time_array, field_array, plot_options)

    # Plot the geometry boundaries
    diagram_tools.plot_geometry_boundaries(coord_array, time_array)
    # Interface target / projectile
    if diagram_tools.data_has_interface():
        diagram_tools.plot_interface(coord_array, time_array)

    # ----------------------------------------------------------
    # Show figure
    # ----------------------------------------------------------
    show_figure(FIG, min_field, max_field)


if __name__ == '__main__':
    N_COLORS = 500

    # ----------------------------------------------------------
    # Read user instructions
    # ----------------------------------------------------------
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("-v", "--verbose", action="store_true", help="Increase program verbosity")
    PARSER.add_argument("field", help="the field to be plotted")
    PARSER.add_argument("case", help="the path to the output repository")
    PARSER.add_argument("-gradient", action="store_true",
                        help="Plot the gradient map instead of the map")
    PARSER.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    ARGS = PARSER.parse_args()

    if ARGS.verbose:
        print("Case : " + ARGS.case)
        print("Field : " + ARGS.field)
        print("Gradient : " + str(ARGS.gradient))
        print("~~~~~~~~~~~~~")

    # ----------------------------------------------------------------
    # Figure creation and settings
    # ----------------------------------------------------------------
    FIG = create_figure()
    # ----------------------------------------------------------
    # Run post processing
    # ----------------------------------------------------------
    run()
