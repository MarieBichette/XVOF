#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
""" 
Compare velocities of free surfaces from simulation with experimental data
USER INPUT :
- The case name to be compared with experiment
(field to be compare is VelocityField for the moment)
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib2tikz import save as tikz_save
from xfv.src.figure_manager.time_figure_tools import VelocityField


project_dir = os.path.split(os.path.dirname(os.path.abspath(os.path.curdir)))[0]
base_dir = os.path.split(project_dir)[0]

test = sys.argv[1]

experimental_file = "{:}.TXT".format(test)
experimental_dir = os.path.join(base_dir, "Documents/These/Experimental_data/Ta/")
exp_data = os.path.join(experimental_dir, experimental_file)

fig = plt.figure(VelocityField.colonne_history)
fig.patch.set_facecolor("white")

fig.suptitle("Vitesse de surface libre", fontsize=20, fontweight='bold')
plt.xlabel("Temps", fontsize=18)
plt.ylabel("Vitesse de la surface libre [m/s]", fontsize=18)



exp = np.loadtxt(exp_data)
plt.plot(exp[:, 0], exp[:, 1], color="black", label="experiment", linestyle="--")

sim = np.loadtxt("CFRAC/{:}/vitesse_surface_libre.dat".format(test))
time = sim[:, 0]
velocity = sim[:, 1]
mask = (velocity > 1)
time_0 = time[np.where(mask)[0][0]]
plt.plot((time - time_0) * 1.e+06, velocity, color="blue", marker=".", label="XFV method", linestyle="-")

sim = np.loadtxt("CFRAC/{:}_pimpose/vitesse_surface_libre.dat".format(test))
time = sim[:, 0]
velocity = sim[:, 1]
mask = (velocity > 1)
time_0 = time[np.where(mask)[0][0]]
plt.plot((time - time_0) * 1.e+06, velocity, color="green", marker=".", label="p=0, S=0", linestyle="-")


plt.legend(loc="best")

tikz_save(os.path.join(base_dir, "Documents/These/Matplotlib2Tikz/experience_ecaillage_{:}.tex".format(test)))
plt.show()
