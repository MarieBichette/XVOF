#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-

from math import pi

import matplotlib.pyplot as plt
import numpy as np
from xvof.element.element1denriched import Element1dEnriched
from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.figure_manager.figure_manager import FigureManager
from xvof.mesh.mesh1denriched import Mesh1dEnriched
from xvof.miscellaneous import geometrical_props, material_props
from xvof.miscellaneous import numerical_props, properties
from xvof.node.node1denriched import Node1dEnriched
from xvof.pressurelaw.constantpressure import ConstantPressure
from xvof.pressurelaw.twostepspressure import TwoStepsPressure
from xvof.rupturecriterion.minimumpressure import MinimumPressureCriterion
from xvof.rupturetreatment.enrichelement import EnrichElement


def print_infos_about_enrichment(mesh, titre="", cells=None, nodes=None):
    has_been_enriched = False
    print "<--- {} --->".format(titre)
    if cells is not None:
        print " Nombre de mailles : {:d}".format(len(cells))
        for cell in cells:
            if isinstance(cell, Element1dEnriched):
#                 indice = cell.indice
#                 cell_g = mesh.cells[indice - 1]
#                 cell_d = mesh.cells[indice + 1]
#                 cell_g.infos()
                cell.infos()
#                 cell_d.infos()
                has_been_enriched = True
    if nodes is not None:
        print " Nombre de noeuds : {:d}".format(len(nodes))
        for node in nodes:
            if isinstance(node, Node1dEnriched):
#                 index = node.index
#                 node_g = mesh.nodes[index - 1]
#                 node_d = mesh.nodes[index + 1]
#                 node_g.infos()
                node.infos()
#                 node_d.infos()
                has_been_enriched = True
    if has_been_enriched:
        raw_input("Poursuivre?")
    print "<--- {} --->".format("-"*len(titre))


def print_infos_about(mesh, titre="", cells=None, nodes=None):
    print "<--- {} --->".format(titre)
    if cells is not None:
        for cell in cells:
            cell.infos()
    if nodes is not None:
        for node in nodes:
            node.infos()
    print "<--- {} --->".format("-"*len(titre))
    raw_input("Poursuivre?")


#  =================================================
#  = PARAMETRES DE LA SIMULATION                   =
TempsFinal = 15.0e-06
PasDeTempsInit = 4.0e-09
PressionInit = 100149.28
EnergieInterneInit = 7.689
RhoInit = 8129.
Section = pi * 0.01 ** 2
EquationEtat = MieGruneisen()
# PChargementGauche = ConstantPressure(-3.5e+09)
PChargementGauche = TwoStepsPressure(15e+09, PressionInit, 2.0e-06)
PChargementDroite = ConstantPressure(PressionInit)
# PChargementDroite = ConstantPressure(-3.5e+09)
CritereRupture = MinimumPressureCriterion(-7.0e+09)
TraitementRupture = EnrichElement(0.5)
Longueur = 10.0e-03
NbrElements = 51
ParamPseudoA = 1.5
ParamPseudoB = 0.2
CFL = 0.35

NbrImages = 150  # 3750
#  =================================================

if __name__ == '__main__':
    #
    time = 0.
    step = 0
    dt = PasDeTempsInit
    dt_crit = 2 * dt
    # ---------------------------------------------#
    #         CREATION DES PROPRIETES              #
    # ---------------------------------------------#
    num_props = numerical_props(ParamPseudoA, ParamPseudoB, CFL)
    mat_props = material_props(PressionInit, EnergieInterneInit, RhoInit, EquationEtat)
    geom_props = geometrical_props(Section)
    props = properties(num_props, mat_props, geom_props)
    # ---------------------------------------------#
    #         CREATION DU MAILLAGE                 #
    # ---------------------------------------------#
    coord_init = np.linspace(0, Longueur, NbrElements + 1)
    vit_init = np.zeros(NbrElements + 1)
    my_mesh = Mesh1dEnriched(props, initial_coordinates=coord_init,
                     initial_velocities=vit_init)
    # ---------------------------------------------#
    #  MISE EN PLACE DU GESTIONNAIRE DE FIGURES    #
    # ---------------------------------------------#
    if (NbrImages != 0):
        delta_t_images = TempsFinal / NbrImages
        my_fig_manager = FigureManager(my_mesh, dump=True, show=True)
        my_fig_manager.populate_figs()
    else:
        delta_t_images = TempsFinal * 2.0
    t_next_image = delta_t_images
    # ---------------------------------------------#
    #         CALCUL DES MASSES NODALES            #
    # ---------------------------------------------#
    print "Calcul de la masse des noeuds :"
    my_mesh.calculer_taille_des_elements()
    my_mesh.calculer_masse_des_noeuds()
    print "=> OK"
    print "LANCEMENT DU CALCUL!"
    while (time < TempsFinal):
#         if(dt_crit < dt):
#             raise SystemExit("Le pas de temps critique est plus petit que le pas de temps")
        print "Itération N°{:<4d} -- Calcul du temps {:15.9g} secondes avec un pas de temps de {:15.9g} secondes"\
            .format(step, time, dt)
#         print "   <- Pas de temps critique = {:15.9g} ->".format(dt_crit)
#         print_infos_about_enrichment(my_mesh, titre="DEBUT DE CYCLE", cells=my_mesh.cells, nodes=my_mesh.nodes)
#         print_infos_about(my_mesh, titre="DEBUT DE CYCLE", cells=my_mesh.cells, nodes=my_mesh.nodes)
        # ---------------------------------------------#
        #         CALCUL DES VITESSES NODALES          #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_vit_noeuds(dt)
#         print_infos_about_enrichment(my_mesh, titre="VITESSES NODALES", nodes=my_mesh.nodes)
        # ---------------------------------------------#
        #         CALCUL DES COORDONNEES NODALES       #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_coord_noeuds(dt)
#         print_infos_about_enrichment(my_mesh, titre="COORDONNEES NODALES", nodes=my_mesh.nodes)
        # ---------------------------------------------#
        #         CALCUL DES VOLUMES DES MAILLES       #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_taille_des_elements(dt)
#         print_infos_about_enrichment(my_mesh, titre="VOLUMES DES MAILLES", cells=my_mesh.cells)
        # ---------------------------------------------#
        #         CALCUL DES DENSITES DES MAILLES      #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_densite_des_elements()
#         print_infos_about_enrichment(my_mesh, titre="DENSITE DES MAILLES", cells=my_mesh.cells)
        # ---------------------------------------------#
        #         CALCUL DES PRESSIONS                 #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_pression_des_elements()
#         print_infos_about_enrichment(my_mesh, titre="PRESSION DES MAILLES", cells=my_mesh.cells)
        # ---------------------------------------------#
        #              RUPTURE                         #
        # ---------------------------------------------#
        my_mesh.get_ruptured_cells(CritereRupture)
        my_mesh.apply_rupture_treatment(TraitementRupture)
#         print_infos_about_enrichment(my_mesh, titre="APRES RUPTURE", cells=my_mesh.cells, nodes=my_mesh.nodes)
        # ---------------------------------------------#
        #         CALCUL DES FORCES NODALES            #
        # ---------------------------------------------#
        my_mesh.calculer_nouvo_force_des_noeuds()
#         print_infos_about_enrichment(my_mesh, titre="FORCES NODALES", nodes=my_mesh.nodes)
        # ---------------------------------------------#
        #         APPLICATION DU CHARGEMENT            #
        # ---------------------------------------------#
        my_mesh.appliquer_pression('gauche', PChargementGauche.evaluate(time))
        my_mesh.appliquer_pression('droite', PChargementDroite.evaluate(time))
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
#         print_infos_about_enrichment(my_mesh, titre="INCREMENTATION", cells=my_mesh.cells, nodes=my_mesh.nodes)
#         dt = min([dt, num_props.cfl * dt_crit])
        time += dt
        step += 1
        # ---------------------------------------------#
        #                GESTION DES SORTIES           #
        # ---------------------------------------------#
        if (time > t_next_image):
            print "Affichage des images"
            my_fig_manager.update_figs("t={:5.4g} us".format(time / 1.e-06))
            t_next_image += delta_t_images
            print "=>OK"

    plt.show()
