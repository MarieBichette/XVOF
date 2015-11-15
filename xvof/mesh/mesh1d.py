#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base définissant un maillage 1d
"""
import numpy as np
from xvof.element.element1d import Element1d
from xvof.mesh.topology1d import Topology1D
from xvof.node.node1d import Node1d


class Mesh1d(object):
    """
    Une classe définissant un maillage 1d
    """
    def __init__(self, properties, initial_coordinates=np.linspace(0, 1, 11),
                 initial_velocities=np.zeros(11)):
        if np.shape(initial_coordinates) != np.shape(initial_velocities):
            message = "Les vecteurs initiaux de vitesse et coordonnées"
            message += " n'ont pas la même taille"
            raise ValueError(message)
        if np.shape(initial_coordinates)[1] != 1:
            message = "Il s'agit d'un maillage 1D à initialiser avec des"
            message += " vecteurs à une dimension!"
            raise ValueError(message)
        ####
        # Création des noeuds
        nbr_nodes = np.shape(initial_coordinates)[0]
        self.nodes = Node1d(nbr_nodes, initial_coordinates, initial_velocities)

        ####
        # Création des mailles
        nbr_cells = nbr_nodes - 1
        self.cells = Element1d(nbr_cells, properties)
        ####
        # Création de la topologie
        self.__topologie = Topology1D(nbr_nodes, nbr_cells)
        ####
        self.__ruptured_cells = []

    def calculer_masse_des_noeuds(self):
        """ Calcul de la masse de chaque noeud"""
        vecteur_nb_noeuds_par_element = np.zeros([nbr_cells, ], type=np.int, order='C')
        vecteur_nb_noeuds_par_element[:] = 2
        self.nodes.calculer_masse_wilkins(self.__topologie, self.cells.masse, vecteur_nb_noeuds_par_element)

    def calculer_nouvo_vit_noeuds(self, delta_t):
        """ Calcul de la nouvelle vitesse de chaque noeud à t+dt"""
        self.nodes.calculer_nouvo_vitesse(delta_t)

    def calculer_nouvo_coord_noeuds(self, delta_t):
        """ Calcul des nouvelles coordonnées de chaque noeud à t+dt"""
        self.nodes.calculer_nouvo_coord(delta_t)

    def calculer_taille_des_elements(self):
        '''
        Calcul de la taille des éléments à t
        '''
        self.cells.computeSize(self.__topologie, self.nodes.xt)

    def calculer_nouvo_taille_des_elements(self, delta_t):
        """ Calcul de la nouvelle taille de chaque élément à t+dt"""
        self.cells.computeNewSize(self.__topologie, self.nodes.xtpdt, delta_t)

    def calculer_nouvo_densite_des_elements(self):
        """ Calcul des nouvelles densités de chaque élément à t+dt"""
        self.cells.computeNewDensity()

    def calculer_nouvo_pression_des_elements(self):
        """ Calcul des nouvelles pressions de chaque élément à t+dt"""
        self.cells.computeNewPressure(mask=self.__ruptured_cells)

    def calculer_nouvo_pseudo_des_elements(self, delta_t):
        """ Calcul de la nouvelle pseudo à t+dt"""
        self.cells.computeNewPseudo(delta_t)

    def calculer_nouvo_force_des_noeuds(self):
        """ Calcul des nouvelles forces de chaque noeud à t+dt"""
        self.nodes.calculer_nouvo_force(self.__topologie, self.cells.pressure.new_value, self.cells.pseudo.new_value)

    def incrementer(self):
        """ Passage au pas de temps suivant"""
        self.nodes.incrementer()
        self.cells.incrementVariables()

    def calculer_nouvo_pdt_critique(self):
        """ Calcul du pas de temps critique """
        self.cells.computeNewTimeStep()
        return self.cells.dt.min()

    def appliquer_pression(self, surface, pression):
        """
        Appliquer une pression donnée sur
        les frontieres gauche ou droite
        """
        if surface.lower() not in ("gauche", "droite"):
            raise(ValueError("Sur la surface <gauche> ou <droite> est possible en 1d!"))
        if (surface.lower() == 'gauche'):
            self.nodes.appliquer_pression(0, pression)
#             self.cells.pressure.new_value[0] = pression
        else:
            self.nodes.appliquer_pression(-1, -pression)
#             self.cells.pressure.new_value[-1] = pression

    @property
    def velocity_t_minus_half_field(self):
        """ Champ de vitesse à t-1/2"""
        return self.nodes.umundemi

    @property
    def velocity_t_plus_half_field(self):
        """ Champ de vitesse à t+1/2"""
        return self.nodes.upundemi

    @property
    def coord_t_field(self):
        """ Champ de position à t"""
        return self.nodes.xt

    @property
    def coord_t_plus_dt_field(self):
        """ Champ de position à t+dt"""
        return self.nodes.xtpdt

    @property
    def coord_elements_field(self):
        """
        Champ de position des éléments à t
        (Moyenne des champs de position à t des noeuds)
        """
        return self.cells.getCoordinates(self.cells.number_of_cells, self.__topologie, self.nodes.xt)

    @property
    def force_field(self):
        """ Champ de force nodale"""
        return self.nodes.force

    @property
    def size_t_field(self):
        """ Tailles des éléments à t"""
        return self.cells.size_t

    @property
    def size_t_plus_dt_field(self):
        """ Tailles des éléments à t"""
        return self.cells.size_t_plus_dt

    @property
    def pressure_t_field(self):
        """ Champ de pression à t"""
        self.cells.pressure.current_value

    @property
    def pressure_t_plus_dt_field(self):
        """ Champ de pression à t+dt"""
        return self.cells.pressure.new_value

    @property
    def rho_t_field(self):
        """ Champ de densité à t"""
        return self.cells.density.current_value

    @property
    def rho_t_plus_dt_field(self):
        """ Champ de densité à t+dt"""
        return self.cells.density.new_value

    @property
    def nrj_t_field(self):
        """ Champ d'énergie interne à t"""
        return self.cells.energy.current_value

    @property
    def pseudo_field(self):
        """ Champ de pseudo """
        return self.cells.pseudo.current_value

    def get_ruptured_cells(self, rupture_criterion):
        """ Liste des mailles endommagées"""
        for ielem in xrange(self.cells.number_of_cells):
            if rupture_criterion.checkCriterion(ielem):
                self.__ruptured_cells.append(ielem)

    def apply_rupture_treatment(self, treatment):
        """
        Application du traitement de rupture sur la liste
        de cells passée en arguments
        """
        for ielem in self.__ruptured_cells:
            treatment.applyTreatment(ielem)
