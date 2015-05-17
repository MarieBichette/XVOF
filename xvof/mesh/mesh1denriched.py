#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant un maillage 1d
"""
from configobj import Section

import numpy as np
from xvof.element.element1d import Element1d
from xvof.node.node1d import Node1d

class Mesh1dEnriched(object):
    """
    Une classe définissant un maillage 1d
    """
    def __init__(self, proprietes, initial_coordinates=np.linspace(0, 1, 11),
                 initial_velocities=np.zeros(11)):
        if np.shape(initial_coordinates) != np.shape(initial_velocities):
            message = "Les vecteurs initiaux de vitesse et coordonnées"
            message += " n'ont pas la même taille"
            raise ValueError(message)
        if len(np.shape(initial_coordinates)) != 1:
            message = "Il s'agit d'un maillage 1D à initialiser avec des"
            message += " vecteurs à une dimension!"
            raise ValueError(message)
        self.__nbr_nodes = np.shape(initial_coordinates)[0]
        self.__nbr_cells = self.__nbr_nodes - 1
        self.__nodes = []
        self.__cells = []
        self.__ruptured_cells = []
        # Création des noeuds
        for n in xrange(self.__nbr_nodes):
            poz = initial_coordinates[n]
            vit = initial_velocities[n]
            nod = Node1d(n, poz_init=np.array([poz]),
                         vit_init=np.array([vit]),
                         section=proprietes.geometric.section)
            self.__nodes.append(nod)
        # Création des éléments
        for m in xrange(self.__nbr_cells):
            elem = Element1d(proprietes, m, [self.__nodes[m], self.__nodes[m + 1]])
            self.__cells.append(elem)
            self.__nodes[m].elements_voisins = [elem]
            self.__nodes[m + 1].elements_voisins = [elem]

    @property
    def nodes(self):
        """ Liste des noeuds """
        return self.__nodes

    @property
    def cells(self):
        """ Liste des éléments """
        return self.__cells

    def calculer_masse_des_noeuds(self):
        """ Calcul de la masse de chaque noeud"""
        for noeud in self.nodes:
            noeud.calculer_masse_wilkins()

    def calculer_nouvo_vit_noeuds(self, delta_t):
        """ Calcul de la nouvelle vitesse de chaque noeud à t+dt"""
        for noeud in self.nodes:
            noeud.calculer_nouvo_vitesse(delta_t)

    def calculer_nouvo_coord_noeuds(self, delta_t):
        """ Calcul des nouvelles coordonnées de chaque noeud à t+dt"""
        for noeud in self.nodes:
            noeud.calculer_nouvo_coord(delta_t)

    def calculer_nouvo_taille_des_elements(self, delta_t):
        """ Calcul de la nouvelle taille de chaque élément à t+dt"""
        for cell in self.cells:
            cell.calculer_nouvo_taille(delta_t)

    def calculer_nouvo_densite_des_elements(self):
        """ Calcul des nouvelles densités de chaque élément à t+dt"""
        for cell in self.cells:
            cell.calculer_nouvo_densite()

    def calculer_nouvo_pression_des_elements(self):
        """ Calcul des nouvelles pressions de chaque élément à t+dt"""
        for cell in self.cells:
            if cell not in self.__ruptured_cells:
                cell.calculer_nouvo_pression()

    def calculer_nouvo_pseudo_des_elements(self, delta_t):
        """ Calcul de la nouvelle pseudo à t+dt"""
        for cell in self.cells:
            cell.calculer_nouvo_pseudo(delta_t)

    def calculer_nouvo_force_des_noeuds(self):
        """ Calcul des nouvelles forces de chaque noeud à t+dt"""
        for noeud in self.nodes:
            noeud.calculer_nouvo_force()

    def incrementer(self):
        """ Passage au pas de temps suivant"""
        for noeud in self.nodes:
            noeud.incrementer()
        for cell in self.cells:
            cell.incrementer()

    def calculer_nouvo_pdt_critique(self):
        """ Calcul du pas de temps critique """
        dts = []
        for cell in self.cells:
            cell.calculer_nouvo_dt()
            dts.append(cell.delta_t)
        return min(dts)

    def appliquer_pression(self, surface, pression):
        """
        Appliquer une pression donnée sur 
        les frontieres gauche ou droite
        """
        if surface.lower() not in ("gauche", "droite"):
            raise(ValueError("Sur la surface <gauche> ou <droite> est possible en 1d!"))
        if (surface.lower() == 'gauche'):
            self.__nodes[0].appliquer_pression(pression)
        else:
            self.__nodes[-1].appliquer_pression(-pression)

    @property
    def velocity_t_minus_half_field(self):
        """ Champ de vitesse à t-1/2"""
        return [node.umundemi for node in self.nodes]

    @property
    def velocity_t_plus_half_field(self):
        """ Champ de vitesse à t+1/2"""
        return [node.upundemi for node in self.nodes]

    @property
    def coord_t_field(self):
        """ Champ de position à t"""
        return [node.coordt for node in self.nodes]

    @property
    def coord_t_plus_dt_field(self):
        """ Champ de position à t+dt"""
        return [node.coordtpdt for node in self.nodes]

    @property
    def coord_elements_field(self):
        """
        Champ de position des éléments à t
        (Moyenne des champs de position à t des noeuds)
        """
        return [cell.coord for cell in self.cells]

    @property
    def force_field(self):
        """ Champ de force nodale"""
        return [node.force for node in self.nodes]

    @property
    def size_t_field(self):
        """ Tailles des éléments à t"""
        return [elem.taille_t for elem in self.cells]

    @property
    def size_t_plus_dt_field(self):
        """ Tailles des éléments à t"""
        return [elem.taille_t_plus_dt for elem in self.cells]

    @property
    def pressure_t_field(self):
        """ Champ de pression à t"""
        return [elem.pression_t for elem in self.cells]

    @property
    def pressure_t_plus_dt_field(self):
        """ Champ de pression à t+dt"""
        return [elem.pression_t_plus_dt for elem in self.cells]

    @property
    def rho_t_field(self):
        """ Champ de densité à t"""
        return [elem.rho_t for elem in self.cells]

    @property
    def rho_t_plus_dt_field(self):
        """ Champ de densité à t+dt"""
        return [elem.rho_t_plus_dt for elem in self.cells]

    @property
    def nrj_t_field(self):
        """ Champ d'énergie interne à t"""
        return [elem.nrj_t for elem in self.cells]

    @property
    def pseudo_field(self):
        """ Champ de pseudo """
        return [elem.pseudo for elem in self.cells]

    def get_ruptured_cells(self, rupture_criterion):
        """ Liste des mailles endommagées"""
        for elem in self.cells:
            if rupture_criterion.checkCriterion(elem):
                self.__ruptured_cells.append(elem)

    def apply_rupture_treatment(self, treatment):
        """
        Application du traitement de rupture sur la liste
        de cells passée en arguments
        """
#         print "Mailles rompues : {}".format(self.__ruptured_cells)
        ruptured_cells = self.__ruptured_cells[:]
        for cell in ruptured_cells:
#             print "-->Traitement de la maille {}".format(cell)
            treatment.applyTreatment(cell, MAILLES=self.__cells,
                                     MAILLES_ROMPUES=self.__ruptured_cells,
                                     NOEUDS=self.__nodes)
#             for cell in self.__cells[cell.indice - 2:cell.indice + 2]:
#                 print cell
#             print "-->self.__ruptured_cells = {}".format(self.__ruptured_cells)
#             raw_input()
