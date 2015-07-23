#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant un maillage 1d
"""

import numpy as np
from xvof.element.element1d import Element1d
from xvof.element.element1denriched import Element1dEnriched
from xvof.mesh.topology1d import Topology1D
from xvof.node.node1d import Node1d

from xvof.utilities.profilingperso import timeit_file

class Mesh1dEnriched(object):
    """
    Une classe définissant un maillage 1d
    """
    def __init__(self, properties, initial_coordinates=np.linspace(0, 1, 11),
                 initial_velocities=np.zeros(11)):
        if np.shape(initial_coordinates) != np.shape(initial_velocities):
            message = "Les vecteurs initiaux de vitesse et coordonnées"
            message += " n'ont pas la même taille"
            raise ValueError(message)
        if len(np.shape(initial_coordinates)) != 1:
            message = "Il s'agit d'un maillage 1D à initialiser avec des"
            message += " vecteurs à une dimension!"
            raise ValueError(message)
        self.__ruptured_cells = []
        ####
        # Création des noeuds
        nodes = []
        for n in xrange(np.shape(initial_coordinates)[0]):
            poz = initial_coordinates[n]
            vit = initial_velocities[n]
            nod = Node1d(poz_init=np.array([poz]),
                         vit_init=np.array([vit]),
                         section=properties.geometric.section)
            nodes.append(nod)
        nbr_nodes = len(nodes)
        ####
        # Création des mailles
        cells = []
        for c in xrange(nbr_nodes - 1):
            cells.append(Element1d(properties))
        ####
        # Création de la topologie
        self.__topologie = Topology1D(nodes, cells)

    @property
    def cells(self):
        '''
        Renvoie la liste des mailles
        '''
        return self.__topologie.cells

    @timeit_file('/tmp/timer.txt')
    def calculer_masse_des_noeuds(self):
        """ Calcul de la masse de chaque noeud"""
        for noeud in self.__topologie.nodes:
            neighbours_cells = self.__topologie._getCellsInContactWithNode(noeud)
            noeud.calculer_masse_wilkins(neighbours_cells)

    @timeit_file('/tmp/timer.txt')
    def calculer_nouvo_vit_noeuds(self, delta_t):
        """ Calcul de la nouvelle vitesse de chaque noeud à t+dt"""
        for noeud in self.__topologie.nodes:
            noeud.calculer_nouvo_vitesse(delta_t)

    @timeit_file('/tmp/timer.txt')
    def calculer_nouvo_coord_noeuds(self, delta_t):
        """ Calcul des nouvelles coordonnées de chaque noeud à t+dt"""
        for noeud in self.__topologie.nodes:
            noeud.calculer_nouvo_coord(delta_t)

    @timeit_file('/tmp/timer.txt')
    def calculer_taille_des_elements(self):
        '''
        Calcul de la taille des éléments à t
        '''
        for cell in self.__topologie.cells:
            nodes = self.__topologie._getNodesBelongingToCell(cell)
            cell.computeSize(nodes)

    @timeit_file('/tmp/timer.txt')
    def calculer_nouvo_taille_des_elements(self, delta_t):
        """ Calcul de la nouvelle taille de chaque élément à t+dt"""
        for cell in self.__topologie.cells:
            nodes = self.__topologie._getNodesBelongingToCell(cell)
            cell.computeNewSize(nodes, delta_t)

    @timeit_file('/tmp/timer.txt')
    def calculer_nouvo_densite_des_elements(self):
        """ Calcul des nouvelles densités de chaque élément à t+dt"""
        for cell in self.__topologie.cells:
            cell.computeNewDensity()

    @timeit_file('/tmp/timer.txt')
    def calculer_nouvo_pression_des_elements(self):
        """ Calcul des nouvelles pressions de chaque élément à t+dt"""
        for cell in self.__topologie.cells:
            if cell not in self.__ruptured_cells:
                cell.computeNewPressure()

    @timeit_file('/tmp/timer.txt')
    def calculer_nouvo_pseudo_des_elements(self, delta_t):
        """ Calcul de la nouvelle pseudo à t+dt"""
        for cell in self.__topologie.cells:
            cell.computeNewPseudo(delta_t)

    @timeit_file('/tmp/timer.txt')
    def calculer_nouvo_force_des_noeuds(self):
        """ Calcul des nouvelles forces de chaque noeud à t+dt"""
        for noeud in self.__topologie.nodes:
            neighbours_cells = self.__topologie._getCellsInContactWithNode(noeud)
            noeud.calculer_nouvo_force(neighbours_cells, isrightboundary=self.__topologie._isRightBoundary(noeud),
                                       isleftboundary=self.__topologie._isLeftBoundary(noeud))

    @timeit_file('/tmp/timer.txt')
    def incrementer(self):
        """ Passage au pas de temps suivant"""
        for noeud in self.__topologie.nodes:
            noeud.incrementer()
        for cell in self.__topologie.cells:
            cell.incrementVariables()

    @timeit_file('/tmp/timer.txt')
    def calculer_nouvo_pdt_critique(self):
        """ Calcul du pas de temps critique """
        dts = []
        for cell in self.__topologie.cells:
            cell.computeNewTimeStep()
            dts.append(cell.delta_t)
        return min(dts)

    @timeit_file('/tmp/timer.txt')
    def appliquer_pression(self, surface, pression):
        """
        Appliquer une pression donnée sur
        les frontieres gauche ou droite
        """
        if surface.lower() not in ("gauche", "droite"):
            raise(ValueError("Sur la surface <gauche> ou <droite> est possible en 1d!"))
        if (surface.lower() == 'gauche'):
            self.__topologie.nodes[0].appliquer_pression(pression)
        else:
            self.__topologie.nodes[-1].appliquer_pression(-pression)

    @property
    def velocity_t_minus_half_field(self):
        """ Champ de vitesse à t-1/2"""
        return [node.umundemi for node in self.__topologie.nodes]

    @property
    def velocity_t_plus_half_field(self):
        """ Champ de vitesse à t+1/2"""
        return [node.upundemi for node in self.__topologie.nodes]

    @property
    def coord_t_field(self):
        """ Champ de position à t"""
        return [node.coordt for node in self.__topologie.nodes]

    @property
    def coord_t_plus_dt_field(self):
        """ Champ de position à t+dt"""
        return [node.coordtpdt for node in self.__topologie.nodes]

    @property
    def coord_elements_field(self):
        """
        Champ de position des éléments à t
        (Moyenne des champs de position à t des noeuds)
        """
        res = []
        for elem in self.__topologie.cells:
            nodes = self.__topologie._getNodesBelongingToCell(elem)
            if isinstance(elem, Element1dEnriched):
                res.append(elem.getLeftPartCoordinates(nodes))
                res.append(elem.getRightPartCoordinates(nodes))
            elif isinstance(elem, Element1d):
                res.append(elem.getCoordinates(nodes))
        return res

    @property
    def force_field(self):
        """ Champ de force nodale"""
        return [node.force for node in self.__topologie.nodes]

    @property
    def size_t_field(self):
        """ Tailles des éléments à t"""
        return [elem.taille_t for elem in self.__topologie.cells]

    @property
    def size_t_plus_dt_field(self):
        """ Tailles des éléments à t"""
        return [elem.taille_t_plus_dt for elem in self.__topologie.cells]

    @property
    def pressure_t_field(self):
        """ Champ de pression à t"""
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.pressure.current_left_value)
                res.append(elem.pressure.current_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.pressure.current_value)
        return res

    @property
    def pressure_t_plus_dt_field(self):
        """ Champ de pression à t+dt"""
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.pressure.new_left_value)
                res.append(elem.pressure.new_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.pressure.new_value)
        return res

    @property
    def rho_t_field(self):
        """ Champ de densité à t"""
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.density.current_left_value)
                res.append(elem.density.current_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.density.current_value)
        return res

    @property
    def rho_t_plus_dt_field(self):
        """ Champ de densité à t+dt"""
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.density.new_left_value)
                res.append(elem.density.new_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.density.new_value)
        return res

    @property
    def nrj_t_field(self):
        """ Champ d'énergie interne à t"""
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.energy.current_left_value)
                res.append(elem.energy.current_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.energy.current_value)
        return res

    @property
    def pseudo_field(self):
        """ Champ de pseudo """
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.pseudo.current_left_value)
                res.append(elem.pseudo.current_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.pseudo.current_value)
        return res

    def get_ruptured_cells(self, rupture_criterion):
        """ Liste des mailles endommagées"""
        for elem in self.__topologie.cells:
            if rupture_criterion.checkCriterion(elem):
                self.__ruptured_cells.append(elem)

    def apply_rupture_treatment(self, treatment):
        """
        Application du traitement de rupture sur la liste
        de cells passée en arguments
        """
        ruptured_cells = self.__ruptured_cells[:]
        for cell in ruptured_cells:
            treatment.applyTreatment(cell, TOPOLOGIE=self.__topologie,
                                     MAILLES_ROMPUES=self.__ruptured_cells)
