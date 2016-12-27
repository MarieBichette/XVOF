# -*- coding: iso-8859-1 -*-
"""
Implementing the Element1dEnriched class
"""
import numpy as np
import os

from xvof.cell import OneDimensionCell
from xvof.data.data_container import DataContainer
from xvof.fields.enrichedfield import from_geometry_to_classic_field, from_geometry_to_enrich_field
from xvof.fields.field import Field


class OneDimensionEnrichedCell(OneDimensionCell):
    """
    A collection of 1d enriched elements
    """

    def __init__(self, number_of_elements):
        super(OneDimensionEnrichedCell, self).__init__(number_of_elements)
        #
        self._fields_manager.moveClassicalToEnrichedFields(number_of_elements)
        #
        self._pos_disc = 0.5  #  La rupture est au milieu de l'élément
        self._fields_manager["taille_gauche"] = Field(number_of_elements, self.size_t * self._pos_disc,
                                                      self.size_t_plus_dt * self._pos_disc)
        self._fields_manager["taille_droite"] = Field(number_of_elements, self.size_t * (1. - self._pos_disc),
                                                      self.size_t_plus_dt * (1. - self._pos_disc))
        print self._fields_manager
        self._classical = np.empty(self.number_of_cells, dtype=np.bool, order='C')
        self._classical[:] = True

    @property
    def classical(self):
        """
        :return: a mask where True indicate a classical cell
        """
        return self._classical

    @property
    def enriched(self):
        """
        :return: a mask where True indicates an enrich cell
        """
        return ~self.classical

    @property
    def left_size(self):
        """
        :return: Left sizes of the enriched elements
        """
        return self._fields_manager['taille_gauche']

    @property
    def right_size(self):
        """
        :return: Right sizes of the enriched elements
        """
        return self._fields_manager['taille_droite']

    @property
    def mass(self):
        """
        :return: Mass of the elements
        """
        return self.size_t * DataContainer().geometric.section * self.density.current_value

    @property
    def pressure_field(self):
        """
        :return: pressure field
        :rtype: list
        :todo: (GP) transform to numpy.array
        """
        p = []
        for i in xrange(self.number_of_cells):
            if self._classical[i]:
                p.append(self.pressure.current_value[i])
            elif self.enriched[i]:
                p.append(self.pressure.current_left_value[i])
                p.append(self.pressure.current_right_value[i])
        return p

    @property
    def density_field(self):
        """
        :return: density field
        :rtype: list
        :todo: (GP) transform to numpy.array
        """
        p = []
        for i in xrange(self.number_of_cells):
            if self._classical[i]:
                p.append(self.density.current_value[i])
            elif self.enriched[i]:
                p.append(self.density.current_left_value[i])
                p.append(self.density.current_right_value[i])
        return p

    @property
    def energy_field(self):
        """
        :return: energy field
        :rtype: list
        :todo: (GP) transform to numpy.array
        """
        p = []
        for i in xrange(self.number_of_cells):
            if self._classical[i]:
                p.append(self.energy.current_value[i])
            elif self.enriched[i]:
                p.append(self.energy.current_left_value[i])
                p.append(self.energy.current_right_value[i])
        return p

    @property
    def artificial_viscosity_field(self):
        """
        :return: artificial viscosity field
        :rtype: list
        :todo: (GP) transform to numpy.array
        """
        p = []
        for i in xrange(self.number_of_cells):
            if self._classical[i]:
                p.append(self.pseudo.current_value[i])
            elif self.enriched[i]:
                p.append(self.pseudo.current_left_value[i])
                p.append(self.pseudo.current_right_value[i])
        return p

    def get_left_part_coordinates(self, topology, nodes_coord):
        """
        :return: coordinates of the left part of enriched elements
        :rtype: numpy.array
        """
        connectivity = topology.nodes_belonging_to_cell
        vec_min = min(nodes_coord[connectivity[self.enriched]])
        vec_coord = vec_min + self.left_size[self.enriched] / 2.
        return vec_coord

    def get_right_part_coordinates(self, topology, nodes_coord):
        """
        :return: coordinates of the right part of enriched elements
        :rtype: numpy.array
        """
        connectivity = topology.nodes_belonging_to_cell
        vec_max = max(nodes_coord[connectivity[self.enriched]])
        vec_coord = vec_max - self.right_size[self.enriched] / 2.
        return vec_coord

    def __str__(self):
        message = "<--ENRICHED CELLS COLLECTION-->" + os.linesep
        message += "Classical elements are:"
        message += str(self.classical)
        message += "Enriched elements are:"
        message += str(self.enriched)
        return message

    def print_infos(self):
        """
        Printing informations about Elements
        """
        message = "{}\n".format(self.__class__)
        message += "==> masse volumique à gauche à t = {}\n". \
            format(self.density.current_left_value)
        message += "==> masse volumique à droite à t = {}\n". \
            format(self.density.current_right_value)
        message += "==> masse volumique classique à t = {}\n". \
            format(self.density.classical_part.current_value)
        message += "==> masse volumique enrichie à t = {}\n". \
            format(self.density.enriched_part.current_value)
        message += "==> masse volumique à gauche à t+dt = {}\n". \
            format(self.density.new_left_value)
        message += "==> masse volumique à droite à t+dt = {}\n". \
            format(self.density.new_right_value)
        message += "==> masse volumique classique à t+dt = {}\n". \
            format(self.density.classical_part.new_value)
        message += "==> masse volumique enrichie à t+dt = {}\n". \
            format(self.density.enriched_part.new_value)
        message += "==> taille à gauche à t = {}\n". \
            format(self.left_size.current_value)
        message += "==> taille à droite à t = {}\n". \
            format(self.right_size.current_value)
        message += "==> taille à gauche à t+dt = {}\n". \
            format(self.left_size.new_value)
        message += "==> taille à droite à t+dt = {}\n". \
            format(self.right_size.new_value)
        message += "==> pression à gauche à t = {}\n". \
            format(self.pressure.current_left_value)
        message += "==> pression à droite à t = {}\n". \
            format(self.pressure.current_right_value)
        message += "==> vitesse du son à gauche à t = {}\n". \
            format(self.sound_velocity.current_left_value)
        message += "==> vitesse du son à droite à t = {}\n". \
            format(self.sound_velocity.current_right_value)
        message += "==> vitesse du son à gauche à t+dt = {}\n". \
            format(self.sound_velocity.new_left_value)
        message += "==> vitesse du son à droite à t+dt = {}\n". \
            format(self.sound_velocity.new_right_value)
        message += "==> énergie à gauche à t = {}\n". \
            format(self.energy.current_left_value)
        message += "==> énergie à droite à t = {}\n". \
            format(self.energy.current_right_value)
        message += "==> énergie à gauche à t+dt = {}\n". \
            format(self.energy.new_left_value)
        message += "==> énergie à droite à t+dt = {}\n". \
            format(self.energy.new_right_value)
        message += "==> pseudo à gauche = {}\n". \
            format(self.pseudo.current_left_value)
        message += "==> pseudo à droite = {}\n". \
            format(self.pseudo.current_right_value)
        print message

    def compute_enriched_elements_new_pressure(self):
        if self.enriched.any():
            mask = self.enriched
            try:
                if self._external_library is not None:
                    energy_left_new, pressure_left_new, sound_velocity_left_new = \
                        (np.array(x) for x in self._compute_new_pressure_with_external_lib(
                            self.density.current_left_value[mask], self.density.new_left_value[mask],
                            self.pressure.current_left_value[mask], self.pseudo.current_left_value[mask],
                            self.energy.current_left_value[mask], self.energy.new_left_value[mask],
                            self.pressure.new_left_value[mask], self.sound_velocity.new_left_value[mask]))
                    energy_right_new, pressure_right_new, sound_velocity_right_new = \
                        (np.array(x) for x in self._compute_new_pressure_with_external_lib(
                            self.density.current_right_value[mask], self.density.new_right_value[mask],
                            self.pressure.current_right_value[mask], self.pseudo.current_right_value[mask],
                            self.energy.current_right_value[mask], self.energy.new_right_value[mask],
                            self.pressure.new_right_value[mask], self.sound_velocity.new_right_value[mask]))
                else:
                    shape = self.energy.new_left_value[mask].shape
                    pressure_left_new = np.zeros(shape, dtype=np.float64, order='C')
                    sound_velocity_left_new = np.zeros(shape, dtype=np.float64, order='C')
                    dummy = np.zeros(shape, dtype=np.float64, order='C')
                    my_variables = {'EquationOfState': DataContainer().material.eos,
                                    'OldDensity': self.density.current_left_value[mask],
                                    'NewDensity': self.density.new_left_value[mask],
                                    'Pressure': (self.pressure.current_left_value[mask] +
                                                 2. * self.pseudo.current_left_value[mask]),
                                    'OldEnergy': self.energy.current_left_value[mask],
                                    'NewPressure': pressure_left_new,
                                    'NewSoundSpeed': sound_velocity_left_new,
                                    'NewDpOverDe': dummy}
                    self._function_to_vanish.setVariables(my_variables)
                    energy_left_new = self._solver.computeSolution(self.energy.current_left_value[mask])
                    self._function_to_vanish.eraseVariables()
                    shape = self.energy.new_right_value[mask].shape
                    pressure_right_new = np.zeros(shape, dtype=np.float64, order='C')
                    sound_velocity_right_new = np.zeros(shape, dtype=np.float64, order='C')
                    dummy = np.zeros(shape, dtype=np.float64, order='C')
                    my_variables = {'EquationOfState': DataContainer().material.eos,
                                    'OldDensity': self.density.current_right_value[mask],
                                    'NewDensity': self.density.new_right_value[mask],
                                    'Pressure': (self.pressure.current_right_value[mask] +
                                                 2. * self.pseudo.current_right_value[mask]),
                                    'OldEnergy': self.energy.current_right_value[mask],
                                    'NewPressure': pressure_right_new,
                                    'NewSoundSpeed': sound_velocity_right_new,
                                    'NewDpOverDe': dummy}
                    self._function_to_vanish.setVariables(my_variables)
                    energy_right_new = self._solver.computeSolution(self.energy.current_right_value[mask])
                    self._function_to_vanish.eraseVariables()
            except ValueError as err:
                raise err
            self.pressure.new_value[self.enriched] = from_geometry_to_classic_field(pressure_left_new,
                                                                                    pressure_right_new)
            self.pressure.new_enr_value[self.enriched] = from_geometry_to_enrich_field(pressure_left_new,
                                                                                       pressure_right_new)
            #
            self.energy.new_value[mask] = from_geometry_to_classic_field(energy_left_new, energy_right_new)
            self.energy.new_enr_value[mask] = from_geometry_to_enrich_field(energy_left_new, energy_right_new)
            #
            self.sound_velocity.new_value[mask] = from_geometry_to_classic_field(sound_velocity_left_new,
                                                                                 sound_velocity_right_new)
            self.sound_velocity.new_enr_value[mask] = from_geometry_to_enrich_field(sound_velocity_left_new,
                                                                                    sound_velocity_right_new)

    def compute_enriched_elements_new_part_size(self, time_step, topologie, vecteur_vitesse_enr_noeud,
                                                vecteur_vitesse_noeuds):
        if self.enriched.any():
            # Calcul des tailles des parties gauches des éléments enrichis
            connectivity = topologie.nodes_belonging_to_cell[self.enriched]
            self.left_size.new_value[self.enriched] = (self.left_size.current_value[self.enriched] +
                                                       (0.5 * (vecteur_vitesse_noeuds[connectivity[:, 1]] -
                                                               vecteur_vitesse_enr_noeud[connectivity[:, 0]]) -
                                                        0.5 * (vecteur_vitesse_noeuds[connectivity[:, 0]] -
                                                               vecteur_vitesse_enr_noeud[connectivity[:, 0]]))
                                                       * time_step)
            self.right_size.new_value[self.enriched] = (self.right_size.current_value[self.enriched] +
                                                        (0.5 * (vecteur_vitesse_noeuds[connectivity[:, 1]] -
                                                                vecteur_vitesse_enr_noeud[connectivity[:, 1]]) -
                                                         0.5 * (vecteur_vitesse_noeuds[connectivity[:, 0]] -
                                                                vecteur_vitesse_enr_noeud[connectivity[:, 1]]))
                                                        * time_step)

    def compute_enriched_elements_new_density(self):
        if self.enriched.any():
            densite_gauche_t_plus_dt = (self.density.current_left_value[self.enriched] *
                                        self.left_size.current_value[self.enriched] / self.left_size.new_value[
                                            self.enriched])
            densite_droite_t_plus_dt = (self.density.current_right_value[self.enriched] *
                                        self.right_size.current_value[self.enriched] / self.right_size.new_value[
                                            self.enriched])
            self.density.new_value[self.enriched] = \
                from_geometry_to_classic_field(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)
            self.density.new_enr_value[self.enriched] = \
                from_geometry_to_enrich_field(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)

    def compute_enriched_elements_new_pseudo(self, delta_t):
        if self.enriched.any():
            rho_t_gauche = self.density.current_left_value[self.enriched]
            rho_t_plus_dt_gauche = self.density.new_left_value[self.enriched]
            cson_t_gauche = self.sound_velocity.current_left_value[self.enriched]
            pseudo_gauche = \
                OneDimensionCell.compute_pseudo(delta_t, rho_t_gauche,
                                                rho_t_plus_dt_gauche,
                                                self.left_size.new_value[self.enriched],
                                                cson_t_gauche,
                                                DataContainer().numeric.a_pseudo, DataContainer().numeric.b_pseudo)

            rho_t_droite = self.density.current_right_value[self.enriched]
            rho_t_plus_dt_droite = self.density.new_right_value[self.enriched]
            cson_t_droite = self.sound_velocity.current_right_value[self.enriched]
            pseudo_droite = \
                OneDimensionCell.compute_pseudo(delta_t, rho_t_droite,
                                                rho_t_plus_dt_droite,
                                                self.right_size.new_value[self.enriched],
                                                cson_t_droite,
                                                DataContainer().numeric.a_pseudo, DataContainer().numeric.b_pseudo)

            self.pseudo.new_value[self.enriched] = \
                from_geometry_to_classic_field(pseudo_gauche, pseudo_droite)
            self.pseudo.new_enr_value[self.enriched] = \
                from_geometry_to_enrich_field(pseudo_gauche, pseudo_droite)

    def compute_enriched_elements_new_time_step(self):
        if self.enriched.any():
            cfl = DataContainer().numeric.cfl
            cfl_pseudo = DataContainer().numeric.cfl_pseudo
            rho_t_gauche = self.density.current_left_value[self.enriched]
            rho_t_plus_dt_gauche = self.density.new_left_value[self.enriched]
            cson_t_plus_dt_gauche = self.sound_velocity.new_left_value[self.enriched]
            pseudo_t_gauche = self.pseudo.current_left_value[self.enriched]
            pseudo_t_plus_dt_gauche = self.pseudo.new_left_value[self.enriched]
            dt_g = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, rho_t_gauche, rho_t_plus_dt_gauche,
                                                      self.left_size.new_value[self.enriched], cson_t_plus_dt_gauche,
                                                      pseudo_t_gauche, pseudo_t_plus_dt_gauche)

            rho_t_droite = self.density.current_right_value[self.enriched]
            rho_t_plus_dt_droite = self.density.new_right_value[self.enriched]
            cson_t_plus_dt_droite = self.sound_velocity.new_right_value[self.enriched]
            pseudo_t_droite = self.pseudo.current_right_value[self.enriched]
            pseudo_t_plus_dt_droite = self.pseudo.new_right_value[self.enriched]
            dt_d = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, rho_t_droite, rho_t_plus_dt_droite,
                                                      self.right_size.new_value[self.enriched], cson_t_plus_dt_droite,
                                                      pseudo_t_droite, pseudo_t_plus_dt_droite)

            self._dt[self.enriched] = dt_g + dt_d  # Bizarre --> A vérifier
