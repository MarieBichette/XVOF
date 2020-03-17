#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
"""
Compare velocities of free surfaces from experimental data
"""

import os
import sys
import matplotlib.pyplot as plt
from xvof.utilities.experimental_data_exploit import ExperimentalData

msg = "Petit programme tout mignon pour tracer la vitesse de surface libre des données expérimentales \n"
msg += "Ce script attend comme  arguments  : \n"
msg += "- une type de matériau (Ta|Cu) \n"
msg += "- -h ou --help pour afficher l'aide"

if len(sys.argv) != 2:
    print msg
    exit(0)

if sys.argv[1] in ["-h", "--help"]:
    print msg
    exit(0)

if sys.argv[1] not in ["Cu", "Ta"]:
    print msg
    exit(0)

experimental_dir = "//home/marie/Documents/These/Experimental_data/{:}/".format(sys.argv[1])
i = 0
color = ["blue", "red", "green", "orange", "purple", "skyblue"]

for file in os.listdir(experimental_dir):
    if not file.startswith("liste"):
        exp_data = os.path.join(experimental_dir, file)
        exp_result = ExperimentalData(exp_data)
        us1, us2 = exp_result.compute_us1_us2()
        print "Calcul us1 = " + us1
        print "Calcul us2 = " + us2
        exp_result.plot_experimental_results(1, color=color[i])
        i += 1
plt.legend()
plt.show()