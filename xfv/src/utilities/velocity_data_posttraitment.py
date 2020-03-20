#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module pour calculer la divergence de la vitesse ou bien reconstruire le champ de vitesse sur les bords de la
discontinuité
"""
import numpy as np


def compute_velocity_divergence_at_time(my_hd, time, cells_id=None):
    """
    Compute the divergence of velocity for the data base my_hd (format hdf5)
    :param my_hd: OutputDataBaseExploit
    :param time: temps de l'extraction de données
    :param cells_ids : index des noeuds si on ne veut sortir que quelques résultats
    :return:
    """
    cell_status = my_hd.extract_field_at_time("CellStatus", time)[:]
    nb_cells = len(cell_status)
    bool_enriched = np.sum(cell_status)  # la somme sera non nulle si au moins une des cells est enrichie

    vitesse = my_hd.extract_field_at_time("ClassicalNodeVelocity", time)[:]
    div_u = np.zeros([nb_cells], order='C')

    # Compute divergence u for classical case
    for i_cell in range(nb_cells):
        node_g_index = int(i_cell)
        node_d_index = int(i_cell) + 1  # marche bien en 1D car topologie facilite l'acquisition des noeuds connectés à la cell
        div_u[i_cell] = (vitesse[node_d_index] - vitesse[node_g_index])

    if bool_enriched:
        # Si mailles rompues : il faut reconstruire la vitesse et coordonnees sur les bords de la disc. et
        # insérer la div sur la partie droite et remplacer le résultat sur la partie gauche
        add_dof_velocity = my_hd.extract_field_at_time("AdditionalNodeVelocity", time)
        epsilon = add_dof_velocity.attrs["discontinuity_position"]

        # tri du tableau des vitesses enrichies pour pouvoir insérer facilement dans le vecteur div_u
        add_dof_velocity = np.sort(add_dof_velocity, 0)
        enriched_cell_id = add_dof_velocity[:, 0]
        add_u_2g = add_dof_velocity[:, 1]
        add_u_1d = add_dof_velocity[:, 2]

        i_stag = 0

        for i_cell_enriched in enriched_cell_id:
            i_cell_enriched = int(i_cell_enriched)
            node_g_enriched = int(i_cell_enriched)
            node_d_enriched = int(i_cell_enriched) + 1
            vitesse_pour_partie_gauche = vitesse[node_g_enriched] * (1. - epsilon) + add_u_2g[0] * epsilon
            vitesse_pour_partie_droite = vitesse[node_d_enriched] * epsilon + add_u_1d[0] * (1. - epsilon)
            div_u[i_cell_enriched + i_stag] = vitesse_pour_partie_gauche - vitesse[node_g_enriched]
            div_u = np.insert(div_u, i_cell_enriched + i_stag + 1, [vitesse[node_d_enriched] - vitesse_pour_partie_droite])
            i_stag += 1

    cell_size = my_hd.get_cells_true_size_at_time(time)
    div_u = div_u / cell_size

    if cells_id is not None:
        div_u = div_u[cells_id]

    return div_u


def compute_velocity_on_discontinuity_borders_at_time(my_hd, time, node_ids=None):
    """
    Compute the velocity field on discontinuity borders for the data base my_hd (format hdf5)
    :param my_hd: OutputDataBaseExploit
    :param case_is_enriched: bool
    :param time : temps de sortie des données
    :param nodes_ids : index des noeuds si on ne veut sortir que quelques résultats
    :return:
    """
    cell_status = my_hd.extract_field_at_time("CellStatus", time)[:]
    bool_enriched = np.sum(cell_status)  # la somme sera non nulle si au moins une des cells est enrichie

    vitesse = my_hd.extract_field_at_time("ClassicalNodeVelocity", time)[:]


    if bool_enriched:
        # Si mailles rompues : il faut reconstruire la vitesse et coordonnees sur les bords de la disc. et
        # insérer la div sur la partie droite et remplacer le résultat sur la partie gauche
        add_dof_velocity = my_hd.extract_field_at_time("AdditionalNodeVelocity", time)
        epsilon = add_dof_velocity.attrs["discontinuity_position"]

        # tri du tableau des vitesses enrichies pour pouvoir insérer facilement dans le vecteur div_u
        add_dof_velocity = np.sort(add_dof_velocity, 0)
        enriched_cell_id = add_dof_velocity[:, 0]
        add_u_2g = add_dof_velocity[:, 1]
        add_u_1d = add_dof_velocity[:, 2]

        i_stag = 0

        for i_cell_enriched in enriched_cell_id:
            node_g_enriched = int(i_cell_enriched)
            node_d_enriched = int(i_cell_enriched) + 1
            vitesse_pour_partie_gauche = vitesse[node_g_enriched] * (1. - epsilon) + add_u_2g[0] * epsilon
            vitesse_pour_partie_droite = vitesse[node_d_enriched] * epsilon + add_u_1d[0] * (1. - epsilon)
            vitesse = np.insert(vitesse, node_g_enriched + i_stag + 1, [vitesse_pour_partie_gauche[0],
                                                                    vitesse_pour_partie_droite[0]])
            i_stag += 2

    if node_ids is not None:
        vitesse = vitesse[node_ids]

    return vitesse