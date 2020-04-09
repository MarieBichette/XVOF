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
    plt.ylabel("Time [$\mu s$]", fontsize=18)
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
        # the_title = "Space time {:} diagram".format(ARGS.field)
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
    final_ruptured_cell_id = np.where(final_cell_status)[0]
    final_ruptured_cell_id = np.sort(final_ruptured_cell_id)
    # ----------------------------------------------------------------
    # Get data to plot the diagram
    # ----------------------------------------------------------------
    diagram_tools = SpaceTimeDiagramTools(path_to_db, ARGS.verbose)
    coord_array, time_array, field_array = \
        diagram_tools.build_xyz_map_for_contourf_plot(ARGS.field, final_ruptured_cell_id)
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

    if len(final_ruptured_cell_id) >= 1:
        # Offset of the cracked cell ids to be conservative with the number of items of final arrays
        ruptured_cell_id_after_offset = \
            final_ruptured_cell_id + list(range(0, len(final_ruptured_cell_id)))
        if ARGS.verbose:
            print("List of cracked cells :" + str(final_ruptured_cell_id))
            print("=> List of cracked cells after offset:" + str(ruptured_cell_id_after_offset))
        # First part from 0 to first discontinuity
        first_left_index = final_ruptured_cell_id[0]
        diagram_tools.plot_section_of_color_map(coord_array, time_array, field_array,
                                                plot_options, end=first_left_index)
        offset = 1
        # From one discontinuity to another
        for i_rupture_index in ruptured_cell_id_after_offset[:-1]:
            right_current = i_rupture_index + 1
            left_next = ruptured_cell_id_after_offset[offset]
            diagram_tools.plot_section_of_color_map(coord_array, time_array, field_array,
                                                    plot_options,
                                                    begin=right_current, end=left_next)
            offset += 1
        # From the last discontinuity to the end
        last_right_index = ruptured_cell_id_after_offset[-1] + 1
        diagram_tools.plot_section_of_color_map(coord_array, time_array, field_array,
                                                plot_options, begin=last_right_index)
    else:
        # Simple plot with no discontinuities
        diagram_tools.plot_section_of_color_map(coord_array, time_array, field_array, plot_options)

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
