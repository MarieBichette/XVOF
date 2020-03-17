#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
"""
Script to compare elastic perfect plastic results from A3 and XVOF
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from xvof.utilities.case_definition import CaseManager
from xvof.figure_manager.time_figure_tools import PressureField, DeviatoricStressField, PseudoViscosityField
from xvof.output_figure.profile_tools import get_error_value, initialize_profile_figure,read_hdf5_file, \
                                            A3_list, plot_field_from_txt_file, get_field_from_txt_file
from xvof.output_manager.outputdatabaseexploit import OutputDatabaseExploit

# Fields to be compared ----------------
field_list = [PressureField, PseudoViscosityField, DeviatoricStressField]

file_name = "all_fields.hdf5"
extension = '.{:}'.format(file_name.split('.')[-1])
# extension = '.txt'
print "extension filefor XVOF data : " + str(extension)

# Error parameters ---------------------
precison_limite = 1.e-03
a3_offset = 2.e-03   # correction pour le décalage de position enregitré
# error_type = 'Absolute'
# error_type = 'Relative'
error_type = 'Absolute_adim'

# -----------------------------------------
# Initialize treatment and legend
# -----------------------------------------
plt.clf()
plt.close()

# -----------------------------------------
# Read user instructions and check entries coherence
# -----------------------------------------
# check number of arguments
if len(sys.argv) != 5:
    raise ValueError("Invalid number of arguments. Expected : \n \
    - the XVOF case to be considered \n \
    - the A3 case for compraison : elasto, epp_120_MPa, hydro_sans_pseudo, hydro_avec_pseudo\n \
    - figure type : profile or (time) evolution \n \
    - a list of time to plot profile at if figure type is profile \n \
    - a list of cell_id if the figure type is evolution \n")

# check case name
physic_case = sys.argv[1].split(',')
A3_case = sys.argv[2]

if A3_case not in ['elasto', 'epp_120MPa', 'hydro_sans_pseudo', 'hydro_avec_pseudo']:
    raise ValueError('the case {:s} is not in : elasto, epp_120MPa, hydro_sans_pseudo, '
                     'hydro_avec_pseudo'.format(A3_case))

# check type of analysis
figure_type = sys.argv[3]
if figure_type not in ['profile', 'error_profile']:
    raise ValueError('figure type expected : profile or error_profile')

# check numerical values (cell id or time for analysis)
time = sys.argv[4].split(',')
for t in time:
    t = float(t)
    if t not in [0., 1.e-06, 2.e-06, 4.e-06, 6.e-06, 8.e-06, 10.e-06]:
        raise ValueError("Data not available for time {:}\n \
        Enter a time in [0., 1.e-06, 2.e-06, 4.e-06, 6.e-06, 8.e-06, 10.e-06]".format(t))

# -----------------------------------------
# Run user instruction (analysis)
# -----------------------------------------
for phys_case in physic_case:
    directory = "//home/marie/Documents/These/export_A3/"  # Find directory of A3 data
    subdir = {'elasto':'RESULT_ELASTO/',
              'epp_120MPa': 'RESULT_EPP_120_MPa/',
              'hydro_sans_pseudo': 'RESULT_HYDRO_sans_pseudo',
              'hydro_avec_pseudo': 'RESULT_HYDRO'}

    if 'profile' in figure_type.split('_'):
        # correction du répertoire pour profile et error_profile : doit aller dans 'profile' pour tous les deux
        directory = os.path.join(directory, subdir[A3_case], 'profile')
    else:
        directory = os.path.join(directory, subdir[A3_case], figure_type)

    # Create case (type Case())from case name
    case = CaseManager().find_case(phys_case)

    # -----------------------------------------
    # Case Treatment
    #  -----------------------------------------
    for field in field_list:
        print field.title
        # -----------------------------------------
        # Field profile and error profile common treatment
        # -----------------------------------------
        # Read data :
        if extension == '.hdf5':
            file_name = "all_fields.hdf5"
            path_to_xvof = os.path.join("./",case.directory_name, file_name)
            my_hd = OutputDatabaseExploit(path_to_xvof)
        for t in time:
            t = float(t)
            # ^^^^^^ Traitement de la bande hdf5 :^^^^^^^^
            if extension == '.hdf5':
                coord, field_value = read_hdf5_file(my_hd, field.title, float(t))
            elif extension == '.txt':
                file_name = "profil_{:s}_at_time_{}_us.txt".format(A3_list[field.title], int(t * 1.e+06))
                if field == PseudoViscosityField: # pour tricher dans le nom des champs dans le cas de la pseudo:
                    file_name = "profil_{:s}_at_time_{}_us.txt".format("artificialviscosity", int(t * 1.e+06))
                path_to_file = os.path.join(case.directory_name, file_name)
                coord, field_value = get_field_from_txt_file(path_to_file)

            # ^^^^^^^Lecture fichier données A3 ^^^^^^
            a3_file_name = "profil_{:s}_at_time_{}_us.txt".format(A3_list[field.title], int(t * 1.e+06))
            path_to_file = os.path.join(directory, a3_file_name)
            coord_a3, field_a3 = get_field_from_txt_file(path_to_file)
            # Recalage
            coord_a3 += -coord_a3[0] + coord[0]

            # -----------------------------------------
            # Plot profile figure ---------------
            # -----------------------------------------
            if figure_type == 'profile':
                # Figure XVOF :
                fig = initialize_profile_figure(field.colonne_history, "{:} profile at time {:} $\mu s$".
                                                format(field.title, t * 1.e+06))
                ax = fig.add_subplot(1, 1, 1)
                ax.plot(coord * 1.e+03, field_value, marker='+', label=case.label)
                # Figure A3
                plot_field_from_txt_file(field.colonne_history, coord_a3, field_a3, offset=0, multiplicateur=1.e+03)
                plt.ylabel(field.label)
                plt.legend(loc='best')
            # -----------------------------------------
            # Plot Error profile figure ---------------
            # -----------------------------------------
            if figure_type == 'error_profile':
                coord_a3 += - a3_offset * 1.e-03  # correction pour le décalage de position enregitré
                # Calcul de l'erreur
                err_x, err_y= get_error_value(coord_a3, field_a3, coord, field_value, error_calcul=error_type.lower())
                # Tracé de la figure
                fig = initialize_profile_figure(field.colonne_history * 10, "Error in {:} profile at time {:} $\mu s$".
                                                format(field.title, t * 1.e+06))
                ax = fig.add_subplot(1, 1, 1)
                ax.plot(err_x, err_y, marker='.', label='error A3 - XVOF')
                ax.plot(err_x, np.ones_like(err_x) * precison_limite, '--', color='red', label='precision limite')
                plt.ylabel('{:} error in {:} field'.format(error_type, field.title))
                plt.legend(loc='best')

print "Done"
plt.show()

