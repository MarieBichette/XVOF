#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
""" 
Compute the energy dissipated by all cohesive zone activated. (aire sous la courbe sigma=f(delta)
"""
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from xvof.utilities.case_definition import CaseManager
from xvof.output_manager.outputdatabaseexploit import OutputDatabaseExploit
from xvof.utilities.velocity_data_posttraitment import compute_velocity_divergence_at_time

# -----------------------------------------
# Initialize treatment and legend
# -----------------------------------------
plt.clf()
plt.close()

# Préparation de la figure de sortie
fig = plt.figure(1)
fig.patch.set_facecolor("white")
fig.suptitle("Energy dissipated", fontsize=20, fontweight='bold')
plt.xlabel("Time [$\mu$s]")
plt.ylabel("Energy dissipated per surface unit[J/m2]")

plt.show(block=False)

# -----------------------------------------
# Read user instructions
# -----------------------------------------
msg = "Petit programme tout mignon pour tracer les différentes sources de dissipation d'énergie \n"
msg += "Ce script attend comme  arguments  : \n"
msg += "- une liste des cas à traiter (séparés par une virgule, sans espace, None par défaut) \n"
msg += "- -h ou --help pour afficher l'aide\n"

if len(sys.argv) != 2:
    print msg
    exit(0)

if sys.argv[1] in ["-h", "--help"]:
    print msg
    exit(0)

case_list = sys.argv[1].split(',')
case_list = [CaseManager().find_case(case) for case in case_list]

# -----------------------------------------
# Lecture de la bande hdf5 pour trouver les données intéressantes :
# -----------------------------------------
for case in case_list:
    my_hd = OutputDatabaseExploit(os.path.join(case.directory_name, "all_fields.hdf5"))

    # On récupère le nombre de mailles
    nb_cells = len(my_hd.extract_field_at_time("CellStatus", my_hd.saved_times[-1]))
    disc_collection_final = np.where(my_hd.extract_field_at_time("CellStatus", my_hd.saved_times[-1]))[0]
    nb_disc = len(disc_collection_final)

    # -----------------------------------------------------
    # Calcul de la dissipation entre t-dt et t (incrément):
    # -----------------------------------------------------

    increment_dissipation_plastique_at_t = np.zeros(my_hd.nb_saved_times)
    increment_dissipation_pseudo_at_t = np.zeros(my_hd.nb_saved_times)
    increment_dissipation_czm_at_t = np.zeros(my_hd.nb_saved_times)
    ouverture_vecteur_temps = np.zeros([my_hd.nb_saved_times, nb_disc])
    force_vecteur_temps = np.zeros([my_hd.nb_saved_times, nb_disc])

    # Dissipation modèle cohésif
    for index_temps in xrange(my_hd.nb_saved_times):
        t = my_hd.saved_times[index_temps]
        # cracked_cell_id, ouverture, force, mask_existance = my_hd.extract_fields_for_cohesive_zone_model(t)
        # On récupère au temps t l'ouverture et la force de chaque discontinuité qui existe
        # (mask_existance pour filtrer). Rq : pas grave si on mélange l'ordre des disc. au cours du calcul. Il s'agit
        # juste de relever les couple force/ouverture pour chacune d'elles
        # ouverture_vecteur_temps[index_temps, :] = ouverture[mask_existance]
        # force_vecteur_temps[index_temps, :] = force[mask_existance]






    for index_temps in xrange(my_hd.nb_saved_times):
        t = my_hd.saved_times[index_temps]
        longueur_mailles = my_hd.get_cells_true_size_at_time(t)
        # Incrémente Dissipation pseudo-viscosité
        pseudo = my_hd.extract_true_field_at_time("ArtificialViscosity", t)[:, 1]
        rho = my_hd.extract_true_field_at_time("Density", t)[:, 1]
        div_u = compute_velocity_divergence_at_time(my_hd, t)
        increment_dissipation_pseudo_cell_by_cell = - pseudo / rho * div_u * longueur_mailles
        # Pour temps=t, on somme les contributions de toutes les mailles pour avoir l'incrément de dissipation dans
        # le maillage
        increment_dissipation_pseudo_at_t[index_temps] = np.sum(increment_dissipation_pseudo_cell_by_cell)


        # Incrémente Dissipation plasticité
        plastic_strain_rate = my_hd.extract_true_field_at_time("PlasticStrainRate", t)[:, 1:]
        pressure = my_hd.extract_true_field_at_time("Pressure", t)[:, 1]
        dev_stress = my_hd.extract_true_field_at_time("DeviatoricStress", t)[:, 1:]
        stress_without_pseudo = np.zeros_like(dev_stress)

        increment_dissipation_plastique_cell_by_cell = np.zeros_like(pressure)

        for i in range(0, 3):
            stress_without_pseudo[:, i] = dev_stress[:, i]
            increment_dissipation_plastique_cell_by_cell[:] = stress_without_pseudo[:, i] * plastic_strain_rate[:, i] * longueur_mailles

        increment_dissipation_plastique_at_t[index_temps] = np.sum(increment_dissipation_plastique_cell_by_cell)

        # Calcul de la dissipation totale :
        dissipation_plastique = np.cumsum(increment_dissipation_plastique_at_t)
        dissipation_pseudo = np.cumsum(increment_dissipation_pseudo_at_t)
        dissipation_totale = dissipation_plastique + dissipation_pseudo


    plt.figure(1)
    plt.plot(my_hd.saved_times, dissipation_plastique, '.-', label="Dissipation plastique")
    plt.plot(my_hd.saved_times, dissipation_pseudo, '.-', label="Dissipation pseudo viscosite")
    plt.plot(my_hd.saved_times, dissipation_plastique + dissipation_pseudo, '.-', label="total")
    plt.legend(loc="best")

    plt.figure(2)
    plt.plot(my_hd.saved_times, dissipation_plastique, '.-', label="Dissipation plastique")
    plt.legend(loc="best")

    plt.figure(3)
    plt.plot(my_hd.saved_times, dissipation_pseudo, '.-', label="Dissipation pseudo viscosite")
    plt.legend(loc="best")



plt.show()