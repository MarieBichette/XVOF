#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base d�finissant un maillage 1d
"""

import numpy as np
from xvof.element.element1d import Element1d
from xvof.element.element1denriched import Element1dEnriched
from xvof.mesh.topology1d import Topology1D
from xvof.node.node1d import Node1d

class Mesh1dEnriched(object):
    """
    Une classe d�finissant un maillage 1d
    """
    def __init__(self, properties, initial_coordinates=np.linspace(0, 1, 11),
                 initial_velocities=np.zeros(11)):
        if np.shape(initial_coordinates) != np.shape(initial_velocities):
            message = "Les vecteurs initiaux de vitesse et coordonn�es"
            message += " n'ont pas la m�me taille"
            raise ValueError(message)
        if len(np.shape(initial_coordinates)) != 1:
            message = "Il s'agit d'un maillage 1D � initialiser avec des"
            message += " vecteurs � une dimension!"
            raise ValueError(message)
        self.__ruptured_cells = []
        ####
        # Cr�ation des noeuds
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
        # Cr�ation des mailles
        cells = []
        for _ in xrange(nbr_nodes - 1):
            cells.append(Element1d(properties))
        ####
        # Cr�ation de la topologie
        self.__topologie = Topology1D(nodes, cells)

    def computeNodesMasses(self):
        """ Calcul de la masse de chaque noeud"""
        for noeud in self.__topologie.nodes:
            neighbours_cells = self.__topologie._getCellsInContactWithNode(noeud)
            noeud.calculer_masse_wilkins(neighbours_cells)

    def computeNewNodesVelocities(self, delta_t):
        """ Calcul de la nouvelle vitesse de chaque noeud � t+dt"""
        for noeud in self.__topologie.nodes:
            noeud.calculer_nouvo_vitesse(delta_t)

    def computeNewNodesCoordinates(self, delta_t):
        """ Calcul des nouvelles coordonn�es de chaque noeud � t+dt"""
        for noeud in self.__topologie.nodes:
            noeud.calculer_nouvo_coord(delta_t)

    def computeCellsSizes(self):
        '''
        Calcul de la taille des �l�ments � t
        '''
        for cell in self.__topologie.cells:
            nodes = self.__topologie._getNodesBelongingToCell(cell)
            cell.computeSize(nodes)

    def computeNewCellsSizes(self, delta_t):
        """ Calcul de la nouvelle taille de chaque �l�ment � t+dt"""
        for cell in self.__topologie.cells:
            nodes = self.__topologie._getNodesBelongingToCell(cell)
            cell.computeNewSize(nodes, delta_t)

    def computeNewCellsDensities(self):
        """ Calcul des nouvelles densit�s de chaque �l�ment � t+dt"""
        for cell in self.__topologie.cells:
            cell.computeNewDensity()

    def computeNewCellsPressures(self):
        """ Calcul des nouvelles pressions de chaque �l�ment � t+dt"""
        for cell in self.__topologie.cells:
            if cell not in self.__ruptured_cells:
                cell.computeNewPressure()

    def computeNewCellsPseudoViscosities(self, delta_t):
        """ Calcul de la nouvelle pseudo � t+dt"""
        for cell in self.__topologie.cells:
            cell.computeNewPseudo(delta_t)

    def computeNewNodesForces(self):
        """ Calcul des nouvelles forces de chaque noeud � t+dt"""
        for noeud in self.__topologie.nodes:
            neighbours_cells = self.__topologie._getCellsInContactWithNode(noeud)
            noeud.calculer_nouvo_force(neighbours_cells, isrightboundary=self.__topologie._isRightBoundary(noeud),
                                       isleftboundary=self.__topologie._isLeftBoundary(noeud))

    def increment(self):
        """ Passage au pas de temps suivant"""
        for noeud in self.__topologie.nodes:
            noeud.incrementer()
        for cell in self.__topologie.cells:
            cell.incrementVariables()

    def computeNewTimeStep(self):
        """ Calcul du pas de temps critique """
        min_dt = 1.0e+09
        for cell in self.__topologie.cells:
            cell.computeNewTimeStep()
            if cell.delta_t < min_dt:
                min_dt = cell.delta_t
        return min_dt

    def applyPressure(self, surface, pressure):
        """
        Appliquer une pressure donn�e sur
        les frontieres gauche ou droite
        """
        if surface.lower() not in ("gauche", "droite"):
            raise(ValueError("Sur la surface <gauche> ou <droite> est possible en 1d!"))
        if (surface.lower() == 'gauche'):
            self.__topologie.nodes[0].applyPressure(pressure)
        else:
            self.__topologie.nodes[-1].applyPressure(-pressure)

    @property
    def velocity_t_minus_half_field(self):
        """ Champ de vitesse � t-1/2"""
        return [node.umundemi for node in self.__topologie.nodes]

    @property
    def velocity_field(self):
        """ Champ de vitesse � t+1/2"""
        return [node.upundemi for node in self.__topologie.nodes]

    @property
    def coord_t_field(self):
        """ Champ de position � t"""
        return [node.coordt for node in self.__topologie.nodes]

    @property
    def nodes_coordinates(self):
        """ Champ de position � t+dt"""
        return [node.coordtpdt for node in self.__topologie.nodes]

    @property
    def cells_coordinates(self):
        """
        Champ de position des �l�ments � t
        (Moyenne des champs de position � t des noeuds)
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
        """ Tailles des �l�ments � t"""
        return [elem.taille_t for elem in self.__topologie.cells]

    @property
    def size_t_plus_dt_field(self):
        """ Tailles des �l�ments � t"""
        return [elem.taille_t_plus_dt for elem in self.__topologie.cells]

    @property
    def pressure_field(self):
        """ Champ de pression � t"""
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
        """ Champ de pression � t+dt"""
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.pressure.new_left_value)
                res.append(elem.pressure.new_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.pressure.new_value)
        return res

    @property
    def density_field(self):
        """ Champ de densit� � t"""
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
        """ Champ de densit� � t+dt"""
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.density.new_left_value)
                res.append(elem.density.new_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.density.new_value)
        return res

    @property
    def energy_field(self):
        """ Champ d'�nergie interne � t"""
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.energy.current_left_value)
                res.append(elem.energy.current_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.energy.current_value)
        return res

    @property
    def pseudoviscosity_field(self):
        """ Champ de pseudo """
        res = []
        for elem in self.__topologie.cells:
            if isinstance(elem, Element1dEnriched):
                res.append(elem.pseudo.current_left_value)
                res.append(elem.pseudo.current_right_value)
            elif isinstance(elem, Element1d):
                res.append(elem.pseudo.current_value)
        return res

    def getRupturedCells(self, rupture_criterion):
        """ Liste des mailles endommag�es"""
        for elem in self.__topologie.cells:
            if rupture_criterion.checkCriterion(elem):
                self.__ruptured_cells.append(elem)

    def applyRuptureTreatment(self, treatment):
        """
        Application du traitement de rupture sur la liste
        de cells pass�e en arguments
        """
        ruptured_cells = self.__ruptured_cells[:]
        for cell in ruptured_cells:
            treatment.applyTreatment(cell, TOPOLOGIE=self.__topologie,
                                     MAILLES_ROMPUES=self.__ruptured_cells)
