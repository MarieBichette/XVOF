#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
""" 
Compare the time evolution of nodes and cell for xfem and reference simulation
"""

import matplotlib.pyplot as plt
import sys
from xvof.figure_manager.time_figure_manager import TimeFigureManager

if len(sys.argv) != 4:
    raise ValueError("Ce script attend comme  arguments le type d'item (node|cell) et le numéro de l'item et le type d'analyse (compare|diff|single)")

item = sys.argv[1]
if not item in ['node', 'cell']:
    raise ValueError("Type d'item {:s} inconnu! Le type doit être soit node soit cell".format(item))
try:
    id_item = int(sys.argv[2])
except ValueError:
    raise ValueError("Le numéro d'item doit être un entier")

analysis = sys.argv[3]
if not analysis in ['single','compare', 'diff']:
    raise ValueError("Type d'analysis {:s} inconnu! Le type doit être soit single soit compare soit diff".format(analysis))

#--------------------------------
simulation = 'XFEM'

fig_manager = TimeFigureManager(simulation , item, id_item )

#---------------------------------
plt.clf()
plt.close()

if analysis == 'single' :
    fig_manager.plot_single_time_figure()
if analysis == 'compare' :
    fig_manager.compare_time_fields()
if analysis == 'diff':
    fig_manager.compare_time_difference_fields()

plt.show()