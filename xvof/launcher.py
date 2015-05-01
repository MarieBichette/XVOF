#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

import numpy as np
from xvof.element.element1d import Element1d
from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.mesh.mesh1d import Mesh1d
from xvof.miscellaneous import geometrical_props, material_props
from xvof.miscellaneous import numerical_props, properties


if __name__ == '__main__':
    print "Création des propriétés"
    equation_detat = MieGruneisen()
    num_props = numerical_props(0.2, 1.0, 0.35)
    mat_props = material_props(100149.28, 7.7, 8129., equation_detat)
    geom_props = geometrical_props(1.0e-06)
    props = properties(num_props, mat_props, geom_props)
    print "=>OK"
    tfinal = 150.0e-06
    time = 0.
    step = 0
    dt = 4.0e-09
    pcharg = 0.5e+09
    coord_init = np.linspace(0, 1.0, 101)
    vit_init = np.zeros(101)
    print "Création du maillage"
    my_mesh = Mesh1d(props, initial_coordinates=coord_init, initial_velocities=vit_init)
    print "=> OK"
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
        # print my_mesh.velocity_t_plus_half_field
        print "=> OK"
        print "Calcul des nouvelles coordonnées des noeuds :"
        my_mesh.calculer_nouvo_coord_noeuds(dt)
        # print my_mesh.coord_t_plus_dt_field
        print "=>OK"        
        print "Calcul des nouvelles tailles des éléments :"
        my_mesh.calculer_nouvo_taille_des_elements()
        # print my_mesh.size_t_plus_dt_field
        print "=> OK"
        print "Calcul des nouvelles densités des éléments :"
        my_mesh.calculer_nouvo_densite_des_elements()
        # print my_mesh.rho_t_plus_dt_field
        print "=> OK"
        print "Calcul des nouvelles pressions des éléments :"
        my_mesh.calculer_nouvo_pression_des_elements()
        # print my_mesh.pressure_t_plus_dt_field
        print "=> OK"
        print "Calcul des nouvelles forces nodales :"
        my_mesh.calculer_nouvo_force_des_noeuds()
        # print "Application du chargement :"
        my_mesh.nodes[0]._force[:] += pcharg * geom_props.section
        my_mesh.nodes[-1]._force[:] += -100149.28 * geom_props.section
        # print my_mesh.force_field
        print "=>OK"
        print "Incrémentation"
        dt_crit = my_mesh.calculer_nouvo_pdt_critique()
        my_mesh.incrementer()
        dt = 0.35 * dt_crit
        time += dt
        step += 1
        # raw_input("Next step?")
