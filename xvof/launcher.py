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
    mat_props = material_props(1.0e+05, 0.0, 8129., equation_detat)
    geom_props = geometrical_props(1.0e-06)
    props = properties(num_props, mat_props, geom_props)
    print "=>OK"
    tfinal = 15.0e-06
    time = 0.
    step = 0.
    dt = 4.0e-09
    pcharg = 0.5e+09
    coord_init = np.linspace(0, 10.0e-03, 101)
    vit_init = np.zeros(101)
    print "Création du maillage"
    my_mesh = Mesh1d(props, initial_coordinates=coord_init, initial_velocities=vit_init)
    print "=> OK"
    print "Node 0 : ", my_mesh.nodes[0], "Elem voisins", \
        map(Element1d.__str__, my_mesh.nodes[0].elements_voisins)
    print "Node 6 : ", my_mesh.nodes[6], "Elem voisins", \
        map(Element1d.__str__, my_mesh.nodes[6].elements_voisins)
    print "Calcul de la masse des noeuds :"
    my_mesh.calculer_masse_des_noeuds()
    print "=> OK"
    print "LANCEMENT DU CALCUL!"
    while (time < tfinal):
        print "Calcul du temps {:15.9g} secondes avec un pas de temps de {:15.9g} secondes"\
            .format(time, dt)
        print "Application du chargement :"
        my_mesh.nodes[0]._force[:] = pcharg * geom_props.section
        print "=>OK"
        print "Calcul de la nouvelle vitesse des noeuds :"
        my_mesh.calculer_nouvo_vit_noeuds(dt)
        print "Vitesse du noeud 0 = ", my_mesh.nodes[0].upundemi
        print "Vitesse du noeud 1 = ", my_mesh.nodes[1].upundemi
        print "=> OK"
        print "Calcul des nouvelles coordonnées des noeuds :"
        my_mesh.calculer_nouvo_coord_noeuds(dt)
        print "Nouvelle coordonnée du noeud 0 =", my_mesh.nodes[0].coordtpdt
        print "Nouvelle coordonnée du noeud 1 =", my_mesh.nodes[1].coordtpdt
        print "Calcul des nouvelles tailes des éléments :"
        my_mesh.calculer_nouvo_taille_des_elements()
        print "Taille de l'élément 0 =", my_mesh.cells[0].taille_t
        print "Nouvelle taille de l'élément 0 =", my_mesh.cells[0].taille_t_plus_dt
        print "=> OK"
        print "Calcul des nouvelles densités des éléments :"
        my_mesh.calculer_nouvo_densite_des_elements()
        print "=> OK"
        print "Calcul des nouvelles pressions des éléments :"
        my_mesh.calculer_nouvo_pression_des_elements()
        print "=> OK"
        print "Calcul des nouvelles forces nodales :"
        my_mesh.calculer_nouvo_force_des_noeuds()
        print "=>OK"
        print "Incrémentation"
        my_mesh.incrementer()
        time += dt
        step += 1