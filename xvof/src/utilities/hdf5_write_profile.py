#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-
"""
Posttreat the hdf5 band to write the true fields for cells and nodes selected
Call with arg :
item : cell/node
id_item : int (list of int possible , entre guillemet et sans espace)
case (optional) : list
"""

import os
import sys

import xvof.src.figure_manager.time_figure_tools as fig_tools
import xvof.src.utilities.case_definition as case_path
from xvof.src.output_figure.hdf5_posttreatment_tools import write_profile_from_db, write_evolution_from_db

# -----------------------------------------
# Read user instructions
# -----------------------------------------
if not len(sys.argv) == 5:
    raise ValueError("Ce script attend comme  arguments  : "
                     "\n les champs à tracer (ex: PressureField, VelocityField, DensityField)"
                     "\n et une liste des cas à traiter"
                     " le type de donnees à écrire : profile ou evolution"
                     " une liste de temps à écrire pour 'profile' ou de item_id pour 'evolution'")

# Field to be analyzed
field_list = sys.argv[1].split(',')

# Cases of interest
case_list = sys.argv[2].split(',')
cases = [getattr(case_path, c) for c in case_list]

analysis = sys.argv[3]

if analysis == 'profile':
    time = sys.argv[4].split(',')
    for t in time:
        t = float(t)
elif analysis == 'evolution':
    item_id = sys.argv[4].split(',')
    item_id = map(int, item_id)
else:
    raise ValueError("""Only evolution or profile are allowed as type of analysis""")
# -----------------------------------------
# Read and save databases for selected cases
# -----------------------------------------
for case in cases:

    # Find database location:
    directory = "./{:}".format(case.directory_name)
    path_to_hdf5_db = os.path.join(directory, "all_fields_bcp_sorties.hdf5")

    if analysis == 'profile':
        for field in field_list:
            field = getattr(fig_tools, field)
            for t in time:
                output_file = "profil_{:}_at_time_{:}_us.txt".format(field.title.lower(), str(int(t*1e06)))
                write_profile_from_db(path_to_hdf5_db, output_file, t, field)
        print "Done for case {:}".format(case.case_name)

    if analysis == 'evolution':
        for field in field_list:
            field = getattr(fig_tools, field)
            for id in item_id:
                output_file = "evolution_at_{:}_{:}.txt".format(field.application, id)
                write_evolution_from_db(path_to_hdf5_db, output_file, field.application, id, field)



