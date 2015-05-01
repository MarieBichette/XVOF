#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-

import matplotlib.animation as anim
import matplotlib.pyplot as plt
import numpy as np
from xvof.element.element1d import Element1d
from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.mesh.mesh1d import Mesh1d
from xvof.miscellaneous import geometrical_props, material_props
from xvof.miscellaneous import numerical_props, properties


TFinal = 15.0e-06
DT_init = 4.0e-09
PChargement = 0.5e+09
Longueur = 0.1
NbrElements = 300
CFL = 0.35

NbrImages = 250


if __name__ == '__main__':
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
    #
    pcharg = PChargement
    #
    coord_init = np.linspace(0, Longueur, NbrElements + 1)
    vit_init = np.zeros(NbrElements + 1)
    my_mesh = Mesh1d(props, initial_coordinates=coord_init, initial_velocities=vit_init)
    #
    delta_t_images = tfinal / NbrImages
    t_next_image = delta_t_images
    print "=> OK"
    print "Initialisation des figures"
    fig_pression = plt.figure()
    ax_pression = fig_pression.add_subplot(111)
    X = np.array(my_mesh.coord_elements_field)
    Y = np.array([1.0e+09 for cell in my_mesh.cells])
    pressure_line, = ax_pression.plot(X, Y, '-')
    ax_pression.set_ylim([0, 1.25 * PChargement])
    ax_pression.set_xlabel('Position [m]')
    ax_pression.set_ylabel('Pression [Pa]')
    ax_pression.set_title('Champ de pression')
    plt.show(block=False)
    print "=>OK"
    print "Calcul de la masse des noeuds :"
    my_mesh.calculer_masse_des_noeuds()
    print "=> OK"
    print "LANCEMENT DU CALCUL!"
    dt_crit = 1e+07
    while (time < tfinal):
        if(dt_crit < dt):
            raise SystemExit("Le pas de temps critique est plus petit que le pas de temps")
        print "Itération N°{:<4d} -- Calcul du temps {:15.9g} secondes avec un pas de temps de {:15.9g} secondes"\
            .format(step, time, dt)
        print "   <- Pas de temps critique = {:15.9g} ->".format(dt_crit)
        print "Calcul de la nouvelle vitesse des noeuds :"
        my_mesh.calculer_nouvo_vit_noeuds(dt)
        print "=> OK"
        print "Calcul des nouvelles coordonnées des noeuds :"
        my_mesh.calculer_nouvo_coord_noeuds(dt)
        print "=>OK"
        print "Calcul des nouvelles tailles des éléments :"
        my_mesh.calculer_nouvo_taille_des_elements()
        print "=> OK"
        print "Calcul des nouvelles densités des éléments :"
        my_mesh.calculer_nouvo_densite_des_elements()
        print "=> OK"
        print "Calcul des nouvelles pressions des éléments :"
        my_mesh.calculer_nouvo_pression_des_elements()
        print "=> OK"
        print "Calcul des nouvelles forces nodales :"
        my_mesh.calculer_nouvo_force_des_noeuds()
        my_mesh.nodes[0]._force[:] += pcharg * geom_props.section
        my_mesh.nodes[-1]._force[:] += -100149.28 * geom_props.section
        print "=>OK"
        if (time > t_next_image):
            print "Affichage des images"
            X = np.array(my_mesh.coord_elements_field)
            Y = np.array(my_mesh.pressure_t_plus_dt_field)
            pressure_line.set_xdata(X)
            pressure_line.set_ydata(Y)
            fig_pression.canvas.draw()
            plt.show(block=False)
            t_next_image += delta_t_images
            print "=>OK"
        print "Incrémentation"
        dt_crit = my_mesh.calculer_nouvo_pdt_critique()
        my_mesh.incrementer()
        dt = CFL * dt_crit
        time += dt
        step += 1
