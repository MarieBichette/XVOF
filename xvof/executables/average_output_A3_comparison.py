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
field = PressureField

# Error parameters ---------------------
precison_limite = 1.e-03

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
    - the first XVOF case to be considered \n \
    - the second XVOF case to be considered \n  "
                     "These two cases will be averaged \n \
    - the A3 case for compraison : elasto, epp_120_MPa, hydro_sans_pseudo, hydro_avec_pseudo\n \
    - a time to plot profile at if figure type is profile")

# check case name
physic_case_one = sys.argv[1]
physic_case_two = sys.argv[2]
A3_case = sys.argv[3]

if A3_case not in ['elasto', 'epp_120MPa', 'hydro_sans_pseudo', 'hydro_avec_pseudo']:
    raise ValueError('the case {:s} is not in : elasto, epp_120MPa, hydro_sans_pseudo, '
                     'hydro_avec_pseudo'.format(A3_case))

# check numerical values (cell id or time for analysis)
time = sys.argv[4]
t = float(time)
if t not in [0., 1.e-06, 2.e-06, 4.e-06, 6.e-06, 8.e-06, 10.e-06]:
    raise ValueError("Data not available for time {:}\n \
    Enter a time in [0., 1.e-06, 2.e-06, 4.e-06, 6.e-06, 8.e-06, 10.e-06]".format(t))

# -----------------------------------------
# Run user instruction (analysis)
# -----------------------------------------

# Traitement cas 1
case_one = CaseManager().find_case(physic_case_one)
file_name = "all_fields.hdf5"
path_to_xvof = os.path.join("./",case_one.directory_name, file_name)
my_hd = OutputDatabaseExploit(path_to_xvof)
coord_1, field_value_1 = read_hdf5_file(my_hd, field.title, t)

# Traitement cas 2
case_two = CaseManager().find_case(physic_case_two)
file_name = "all_fields.hdf5"
path_to_xvof = os.path.join("./",case_two.directory_name, file_name)
my_hd = OutputDatabaseExploit(path_to_xvof)
coord_2, field_value_2 = read_hdf5_file(my_hd, field.title, t)

# Moyenne des deux cas XVOF :
averaged_coord = 0.5 * (coord_1 + coord_2)
averaged_field_value = 0.5 * (field_value_1 + field_value_2)

# Cas A3
directory = "//home/marie/Documents/These/export_A3/"  # Find directory of A3 data
subdir = {'elasto':'RESULT_ELASTO/',
          'epp_120MPa': 'RESULT_EPP_120_MPa/',
          'hydro_sans_pseudo': 'RESULT_HYDRO_sans_pseudo',
          'hydro_avec_pseudo': 'RESULT_HYDRO'}
directory = os.path.join(directory, subdir[A3_case], 'profile')
a3_file_name = "profil_{:s}_at_time_{}_us.txt".format(A3_list[field.title], int(t * 1.e+06))
path_to_file = os.path.join(directory, a3_file_name)
coord_a3, field_a3 = get_field_from_txt_file(path_to_file)
# Recalage
coord_a3 += -coord_a3[0] + averaged_coord[0]

# -----------------------------------------
# Plot profile figure ---------------
# -----------------------------------------
# Figure XVOF :
fig = initialize_profile_figure(field.colonne_history, "{:} profile at time {:} $\mu s$".
                                format(field.title, t * 1.e+06))
ax = fig.add_subplot(1, 1, 1)
ax.plot(averaged_coord * 1.e+03, averaged_field_value, marker='+', label="average")
ax.plot(coord_1 * 1.e+03, field_value_1, marker='+', label=case_one.label)
ax.plot(coord_2 * 1.e+03, field_value_2, marker='+', label=case_two.label)
# Figure A3
plot_field_from_txt_file(field.colonne_history, coord_a3, field_a3, offset=0, multiplicateur=1.e+03)
plt.ylabel(field.label)
plt.legend(loc='best')

print "Done"
plt.show()

