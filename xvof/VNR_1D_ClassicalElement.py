#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-

import matplotlib.pyplot as plt
import numpy as np

from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.figure_manager.figure_manager import FigureManager

from xvof.mesh.mesh1d import Mesh1d
from xvof.miscellaneous import geometrical_props, material_props
from xvof.miscellaneous import numerical_props, properties

#  =================================================
#  = PARAMETRES DE LA SIMULATION                   =
TFinal = 15.0e-06
DT_init = 4.0e-09
PChargement = 5.0e+09
Longueur = 0.1
NbrElements = 300
CFL = 0.35

NbrImages = 25
#  =================================================

if __name__ == '__main__':
    # ---------------------------------------------#
    #         CREATION DE L'ENVIRONNEMENT          #
    # ---------------------------------------------#
    print "Création des propriétés"
    equation_detat = MieGruneisen()
    num_props = numerical_props(0.2, 1.0, 0.35)
    mat_props = material_props(100149.28, 7.7, 8129., equation_detat)
    geom_props = geometrical_props(1.0e-06)
    props = properties(num_props, mat_props, geom_props)
    print "=>OK"
    print "Création du maillage"
    tfinal = TFinal
    time = 0.
    step = 0
    dt = DT_init
    dt_crit = 2 * dt
    pcharg = PChargement
    coord_init = np.linspace(0, Longueur, NbrElements + 1)
    vit_init = np.zeros(NbrElements + 1)
    my_mesh = Mesh1d(props, initial_coordinates=coord_init,
                     initial_velocities=vit_init)
    my_fig_manager = FigureManager(my_mesh)
    #
    delta_t_images = tfinal / NbrImages
    t_next_image = delta_t_images
    print "=> OK"
    print "Initialisation des figures"
    my_fig_manager.populate_figs()
    print "=>OK"
    # ---------------------------------------------#
    #         CALCUL DES MASSES NODALES            #
    # ---------------------------------------------#
    print "Calcul de la masse des noeuds :"
    my_mesh.calculer_masse_des_noeuds()
    print "=> OK"
    print "LANCEMENT DU CALCUL!"
    while (time < tfinal):
        if(dt_crit < dt):
            raise SystemExit("Le pas de temps critique est plus petit que le pas de temps")
        print "Itération N°{:<4d} -- Calcul du temps {:15.9g} secondes avec un pas de temps de {:15.9g} secondes"\
            .format(step, time, dt)
        print "   <- Pas de temps critique = {:15.9g} ->".format(dt_crit)
        # ---------------------------------------------#
        #         CALCUL DES VITESSES NODALES          #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_vit_noeuds(dt)
        # ---------------------------------------------#
        #         CALCUL DES COORDONNEES NODALES       #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_coord_noeuds(dt)
        # ---------------------------------------------#
        #         CALCUL DES VOLUMES DES MAILLES       #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_taille_des_elements()
        # ---------------------------------------------#
        #         CALCUL DES DENSITES DES MAILLES      #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_densite_des_elements()
        # ---------------------------------------------#
        #         CALCUL DES PRESSIONS                 #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_pression_des_elements()
        # ---------------------------------------------#
        #         CALCUL DES FORCES NODALES            #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_force_des_noeuds()
        # ---------------------------------------------#
        #         APPLICATION DU CHARGEMENT            #
        # ---------------------------------------------#
        my_mesh.nodes[0]._force[:] += pcharg * geom_props.section
        my_mesh.nodes[-1]._force[:] += -100149.28 * geom_props.section
        # ---------------------------------------------#
        #         CALCUL DU PAS DE TEMPS CRITIQUE      #
        # ---------------------------------------------#
        dt_crit = my_mesh.calculer_nouvo_pdt_critique()
        # ---------------------------------------------#
        #         CALCUL DE LA PSEUDOVISCOSITE         #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_pseudo_des_elements(dt)
        # ---------------------------------------------#
        #                INCREMENTATION                #
        # ---------------------------------------------#
        my_mesh.incrementer()
        dt = CFL * dt_crit
        time += dt
        step += 1
        # ---------------------------------------------#
        #                GESTION DES SORTIES           #
        # ---------------------------------------------#
        if (time > t_next_image):
            print "Affichage des images"
            my_fig_manager.update_figs()
            t_next_image += delta_t_images
            print "=>OK"

    plt.show()
