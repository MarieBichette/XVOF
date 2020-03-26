#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Plot the profil evolution of nodes and cell
"""

import argparse
import pathlib
import numpy as np
import matplotlib.pyplot as plt

from xfv.post_processing.tools.hdf5_postprocessing_tools import get_field_profile_at_time
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit

# ----------------------------------------------------------
# Read instructions
# ----------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", action="store_true", help="Increase program verbosity")
parser.add_argument("field", help="the field to be plotted")
parser.add_argument("time", type=float, action='append', nargs='+',
                    help="the list of times to look at")
parser.add_argument("-case", action='append', nargs='+',
                    help="the path to the output repository")
parser.add_argument("--output_filename", default="all_fields.hdf5",
                    help="the name of the output hdf5 band (default = all_fields.hdf5)")
args = parser.parse_args()

if args.case is None:
    raise ValueError("At least one case is needed. Use -case to specify cases to plot")

if args.time is None:
    raise ValueError("At least one time is needed. Use -time to specify times to plot")

if args.verbose:
    print("Cases : ")
    print(args.case)
    print("~~~~~~~~~~~~~")

field = args.field

# ----------------------------------------------------------
# Prepare figure
# ----------------------------------------------------------
plt.figure(1)
plt.title("{:} profile".format(field), fontweight='bold')
plt.xlabel("Position [mm]")
plt.ylabel("Free surface velocity [m/s]")

# ----------------------------------------------------------
# Plot field profile for each case
# ----------------------------------------------------------
for case in args.case[0]:
    path_to_db = pathlib.Path.cwd().joinpath("..", "tests", case, args.output_filename)
    if args.verbose:
        print("Path to database : {:}".format(path_to_db))
        print("Read field " + field + " in database... ")
    hd_band = OutputDatabaseExploit(path_to_db)

    for time in args.time:
        # Read hdf5 band :
        coord, field_value = get_field_profile_at_time(hd_band, field, time)
        if args.verbose:
            print("Done !")
            print("~~~~~~~~~~~~~")
        # Plot field :
        plt.plot(coord * 1.e+03, field_value, label=case)

# ----------------------------------------------------------
# Show figure
# ----------------------------------------------------------
plt.legend(loc="best")
plt.show()
