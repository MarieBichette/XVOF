#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
"""
tools to plot field profile and comparison between profiles
"""

import matplotlib.pyplot as plt
import numpy as np

# Equivalence field name between fields from xvof and fields from A3
A3_list = {'Pressure':'pressure',
           'DeviatoricStress': 'deviatoricstresstensor_1',
           'EquivalentPlasticStrainRate': 'equivalentplasticstrainrate',
           'InternalEnergy': 'internalenergy',
           'ArtificialViscosity': 'scalarpseudo',
           'Density': 'density'}


def plot_field_from_txt_file(fig_id, x, y, offset=0., multiplicateur=1.):
    """
    read the A3 result in appropriate file and plot
    :param fig_id : id where to plot
    :param x : abscisse
    :param y : ordonnée
    :param offset : décalage pour recaler en temps
    :param multiplicateur : multiplie l'Ã©chelle des abscisses
    """
    my_legend = "code A3"
    ax = plt.figure(fig_id).add_subplot(1, 1, 1)
    ax.plot((x - offset) * multiplicateur, y, marker='.', label=my_legend, color='red')


def get_field_from_txt_file(path):
    """
    Read the file where a3 results are stored
    :param path: path to a3 data file
    :return: x, y
    """
    res = np.loadtxt(path, dtype='float')
    x = res[:, 0]
    y = res[:, 1]
    return x, y

def get_error_value(coord_ref, field_ref, coord, field_value, error_calcul='absolute'):
    """
    Compute the error in profile
    :param coord_ref: coordinates (reference)
    :param field_ref: field value (reference)
    :param coord: coordinates from xvof
    :param field_value: field value from xvof
    :param error_calcul : 'absolute' or 'relative'
    :return: err_x, err_y, interpolated coordinates and error in field value
    """
    # Interpolation sur une même grille :
    print "Error computation requires interpolation on a common fixed grid"
    borne_left = max(min(coord_ref), min(coord))  # boundary of the left part of interpolation domain
    borne_right = min(max(coord_ref), max(coord))    # boundary of the right part of interpolation domain
    x_interp = np.linspace(borne_left, borne_right, num=1000)  # create the fix x-grid for interpolation

    y_interp_a3 = np.interp(x_interp, coord_ref, field_ref, left=-1, right=-1)
    y_interp = np.interp(x_interp, coord, field_value, left=-1, right=-1)

    print "Returning {:} error".format(error_calcul)
    # Calcul de l'erreur absolue
    if error_calcul == 'absolute':
        return x_interp, np.abs(y_interp - y_interp_a3)

    if error_calcul == 'relative':
        return x_interp, np.abs(y_interp - y_interp_a3)/np.abs(y_interp_a3)

    if error_calcul == 'absolute_adim':
        return x_interp, np.abs(y_interp - y_interp_a3) / np.max(np.abs(y_interp_a3))


def initialize_profile_figure(id_fig, figure_suptitle):
    fig = plt.figure(id_fig)
    fig.patch.set_facecolor("white")
    fig.suptitle(figure_suptitle, fontsize=20, fontweight='bold')
    plt.xlabel("Position [$mm$]", fontsize=18)
    return fig


def read_hdf5_file(hdf5_band, field_title, t):
    data = hdf5_band.extract_true_field_at_time(field_title, t)
    coord = data[:, 0]
    field_value = data[:, 1]
    return coord, field_value


def plot_field_with_rectangle_shapes(cells_coordinates, cells_size, field_value, ax):
    """
    :param cells_coordinates:
    :param cells_size:
    :param field_value:
    :param ax: subplot to plot the rectangles in
    :return:
    """
    for i in range(0, len(cells_coordinates)-1):
        x_rect = np.array([cells_coordinates[i] - cells_size[i]/2,
                           cells_coordinates[i] - cells_size[i]/2,
                           cells_coordinates[i] + cells_size[i]/2,
                           cells_coordinates[i] + cells_size[i]/2])
        y_rect = np.array([0, field_value[i], field_value[i], 0])
        ax.plot(x_rect * 1.e+03, y_rect, color="red")
        ax.fill_between(x_rect, y_rect, where=(x_rect >= min(x_rect)) & (x_rect <= max(x_rect)), facecolor="orange")