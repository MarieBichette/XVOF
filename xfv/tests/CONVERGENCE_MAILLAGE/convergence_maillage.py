#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
"""
Posttrait data pour faire convergence en maillage
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib2tikz

from xfv.output_manager.outputdatabaseexploit import OutputDatabaseExploit


####################################################
# Definition de l'erreur
####################################################

def compute_error(cell_id, data, ref_data):
    """
    Norme infinie de l'�cart entre data et ref_data � un temps donn�
    :param cell_id:
    :param data:
    :param ref_data:
    :return:
    """

    # if len(data) != len(ref_data):
    #     raise ValueError("""Impossible de calculer la norme de l'erreur si data et ref_data n'ont pas la m�me taillle""")

    global_error = 0.

    size_ref = len(ref_data)
    size_data = len(data)

    factor = size_ref / size_data

    ref_max = np.max(np.abs(ref_data))

    for cell_id in range(size_data):
        x_target = data[cell_id * factor]
        x_ref = ref_data[cell_id]

        if abs(x_ref) > 0.01 * ref_max:  # pour �viter les cas pathologiques o� le champ = 0 (=> diviser par 0)
            error = abs((x_target - x_ref) / x_ref)
            global_error += error

    return global_error


####################################################
# Definition du temps d'observation
####################################################
def find_time_min_density(ref_hd):
    """
    Identification du temps o� le choc arrive sur la derni�re maille de la ref = equiv de la maille rompue
    :param ref_hd:
    :return:
    """
    min_time = -1
    global_pressure_min = 1.e+10
    for t in ref_hd.saved_times:
        pressure_at_t = ref_hd.extract_true_field_at_time("Pressure", t)[:, 1]
        p_min = np.min(pressure_at_t)
        if p_min < global_pressure_min:
            min_time = t
            global_pressure_min = p_min
    return min_time


####################################################
# Initialisation des constantes
####################################################
field = "Pressure"

####################################################
# Initialisation des tableaux contenant les r�sultats
####################################################
maillage = np.arange(4)
print("maillage = ")
print(maillage)

taille_maillage = np.array([1.e-6, 3.e-6, 9.e-6, 27.e-6, 81.e-6])

erreur_maillage_ref = np.zeros_like(taille_maillage)
erreur_maillage_somme = np.zeros_like(taille_maillage)
erreur_maillage_menouillard = np.zeros_like(taille_maillage)


####################################################
# Lecture des fichiers hdf5 pour r�f�rence
####################################################
dir_to_db = "Ref_maillage_0"
path_to_db = os.path.join(dir_to_db, "all_fields.hdf5")
ref_hd = OutputDatabaseExploit(path_to_db)
# time = find_time_min_density(ref_hd)
time = 0.71e-6  # environ dans le choc...
ref_data = ref_hd.extract_true_field_at_time(field, time)[:, 1]

####################################################
# Lecture des fichiers hdf5
####################################################
for i in maillage:
    i += 1

########################################################################

    print("Maillage de r�f�rence : ")
    dir_to_db = "Ref_maillage_" + str(i)
    path_to_db = os.path.join(dir_to_db, "all_fields.hdf5")
    my_hd = OutputDatabaseExploit(path_to_db)
    data = my_hd.extract_true_field_at_time(field, time)[:, 1]

    ####################################################
    # Filtrage ? Interpolation ?
    ####################################################

    ####################################################
    # Calcul de l'erreur
    ####################################################
    erreur_maillage_ref[i] = compute_error(-1, data, ref_data)

########################################################################

    print("Maillage de somme : ")
    dir_to_db = "Enr_somme_maillage_" + str(i)
    path_to_db = os.path.join(dir_to_db, "all_fields.hdf5")
    my_hd = OutputDatabaseExploit(path_to_db)
    data = my_hd.extract_true_field_at_time(field, time)[:, 1]

    ####################################################
    # Filtrage ? Interpolation ?
    ####################################################

    ####################################################
    # Calcul de l'erreur
    ####################################################
    erreur_maillage_somme[i] = compute_error(-1, data, ref_data)

########################################################################
    print("Maillage de menouillard : ")
########################################################################
    dir_to_db = "Enr_menouillard_maillage_" + str(i)
    path_to_db = os.path.join(dir_to_db, "all_fields.hdf5")
    my_hd = OutputDatabaseExploit(path_to_db)
    data = my_hd.extract_true_field_at_time(field, time)[:, 1]

    ####################################################
    # Filtrage ? Interpolation ?
    ####################################################

    ####################################################
    # Calcul de l'erreur
    ####################################################
    erreur_maillage_menouillard[i] = compute_error(-1, data, ref_data)


####################################################
# Prise d'ordre
####################################################
y05 = 1.e-4
y1 = 1.e3
y2 = 1.e5
ordre05 = y05 * np.sqrt(taille_maillage)
ordre1= y1 * taille_maillage
ordre2 = y2 * taille_maillage * taille_maillage

####################################################
# Trac� des donn�es
####################################################
fig = plt.figure(1)
fig.patch.set_color("white")
fig.suptitle("Etude de convergence a t = {:}".format(time), fontsize=20, fontweight="bold")
plt.xlabel("Taille de maille [m]")
plt.ylabel("Erreur relative sur {:}".format(field))

# trac� enr somme
plt.plot(taille_maillage, erreur_maillage_somme, linestyle="-", marker=".", color="red", label="enr somme")

# trac� enr menouillard
plt.plot(taille_maillage, erreur_maillage_menouillard, linestyle="-", marker=".", color="green", label="enr menouillard")

# trac� ref
plt.plot(taille_maillage, erreur_maillage_ref, linestyle="-", marker=".", color="blue", label="ref")

# trac� ordre 0.5
#plt.plot(taille_maillage, ordre05, color="black", label="ordre 0.5")

# trac� ordre 1
plt.plot(taille_maillage, ordre1, color="orange", label="ordre 1")

# trac� ordre 2
plt.plot(taille_maillage, ordre2, color="black", label="ordre 2")

plt.loglog()
plt.grid(which="both")
plt.legend()
plt.show()