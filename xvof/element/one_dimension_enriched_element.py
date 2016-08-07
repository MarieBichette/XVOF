# -*- coding: iso-8859-1 -*-
"""
Implementing the Element1dEnriched class
"""
import ctypes
import numpy as np

from xvof.element import OneDimensionElement
from xvof.fields.enrichedfield import from_geometry_to_classic_field, from_geometry_to_enrich_field
from xvof.fields.field import Field


class OneDimensionEnrichedOneDimensionElement(OneDimensionElement):
    """
    A collection of 1d enriched elements
    """
    def __init__(self, number_of_elements, properties):
        super(OneDimensionEnrichedOneDimensionElement, self).__init__(number_of_elements, properties)
        #
        self._fields_manager.moveClassicalToEnrichedFields(number_of_elements)
        #
        self._pos_disc = 0.5  #  La rupture est au milieu de l'élément
        self._fields_manager["taille_gauche"] = Field(number_of_elements, self.size_t * self._pos_disc,
                                                      self.size_t_plus_dt * self._pos_disc)
        self._fields_manager["taille_droite"] = Field(number_of_elements, self.size_t * (1. - self._pos_disc),
                                                      self.size_t_plus_dt * (1. - self._pos_disc))
        print self._fields_manager
        self._classical = np.empty(self._shape, dtype=np.bool, order='C')
        self._classical[:] = True

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
    def _enriched(self):
        """
        :return: Boolean mask of the enriched elements
        """
        return ~self._classical

    @property
    def mass(self):
        """
        :return: Mass of the elements
        """
        return self.size_t * self.proprietes.geometric.section * self.density.current_value

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
            elif self._enriched[i]:
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
            elif self._enriched[i]:
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
            elif self._enriched[i]:
                p.append(self.energy.current_left_value[i])
                p.append(self.energy.current_right_value[i])
        return p

    @property
    def pseudoviscosity_field(self):
        """
        :return: artificial viscosity field
        :rtype: list
        :todo: (GP) transform to numpy.array
        """
        p = []
        for i in xrange(self.number_of_cells):
            if self._classical[i]:
                p.append(self.pseudo.current_value[i])
            elif self._enriched[i]:
                p.append(self.pseudo.current_left_value[i])
                p.append(self.pseudo.current_right_value[i])
        return p

    def getLeftPartCoordinates(self, topology, nodes_coord):
        """
        :return: coordinates of the left part of enriched elements
        :rtype: numpy.array
        """
        connectivity = np.array(topology._nodes_belonging_to_cell)
        vec_min = min(nodes_coord[connectivity[self._enriched]])
        vec_coord = vec_min + self.left_size[self._enriched] / 2.
        return vec_coord

    def getRightPartCoordinates(self, topology, nodes_coord):
        """
        :return: coordinates of the right part of enriched elements
        :rtype: numpy.array
        """
        connectivity = np.array(topology._nodes_belonging_to_cell)
        vec_max = max(nodes_coord[connectivity[self._enriched]])
        vec_coord = vec_max - self.right_size[self._enriched] / 2.
        return vec_coord

    def __str__(self):
        message = "ENRICHED ELEMENT {:4d} ".format(self._index)
        return message

    def printInfos(self):
        """
        Printing informations about Elements
        """
        message = "{} {:4d}\n".format(self.__class__, self._index)
        message = "==> masse volumique à gauche à t = {}\n". \
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
            format(self._taille_gauche_t)
        message += "==> taille à droite à t = {}\n". \
            format(self._taille_droite_t)
        message += "==> taille à gauche à t+dt = {}\n". \
            format(self._taille_gauche_t_plus_dt)
        message += "==> taille à droite à t+dt = {}\n". \
            format(self._taille_droite_t_plus_dt)
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

    def computeNewPressure(self):
        """
        Calcul du triplet energie, pression, vitesse du son
        au pas de temps suivant
        Formulation v-e
        """
        # -----------------------------
        # Pression éléments classiques
        # -----------------------------
        super(OneDimensionEnrichedOneDimensionElement, self).computeNewPressure(mask=self._classical)
        if self._enriched.any():
            # -----------------------------
            # Pression partie gauche
            # -----------------------------
            density_current_value = self.density.current_left_value[self._enriched]
            density_new_value = self.density.new_left_value[self._enriched]
            pressure_current_value = self.pressure.current_left_value[self._enriched]
            pseudo_current_value = self.pseudo.current_left_value[self._enriched]
            energy_current_value = self.energy.current_left_value[self._enriched]
            energy_new_value = self.energy.new_left_value[self._enriched]
            shape = energy_new_value.shape
            nbr_cells_to_solve = shape[0]
            solution_value = np.zeros(shape, dtype=np.float64, order='C')
            new_pressure_value = np.zeros(shape, dtype=np.float64, order='C')
            new_vson_value = np.zeros(shape, dtype=np.float64, order='C')
            try:
                if self._external_library is not None:
                    pb_size = ctypes.c_int()
                    #
                    old_density = density_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    new_density = density_new_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    tmp = (pressure_current_value + 2. * pseudo_current_value)
                    pressure = tmp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    old_energy = energy_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    pb_size.value = nbr_cells_to_solve
                    solution = solution_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    new_pressure = new_pressure_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    new_vson = new_vson_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    self._computePressureExternal(old_density, new_density, pressure, old_energy, pb_size,
                                                  solution, new_pressure, new_vson)
                    energy_left_new = solution[0:nbr_cells_to_solve]
                    pressure_left_new = new_pressure[0:nbr_cells_to_solve]
                    sound_velocity_left_new = new_vson[0:nbr_cells_to_solve]
                else:
                    my_variables = {'EquationOfState': self.proprietes.material.eos,
                                    'OldDensity': density_current_value,
                                    'NewDensity': density_new_value,
                                    'Pressure': pressure_current_value + 2. * pseudo_current_value,
                                    'OldEnergy': energy_current_value}
                    self._function_to_vanish.setVariables(my_variables)
                    solution = self._solver.computeSolution(energy_current_value)
                    new_pressure_value, _, new_vson_value = \
                        self.proprietes.material.eos.solveVolumeEnergy(1. / density_new_value, solution)
                    energy_left_new = solution
                    pressure_left_new = new_pressure_value
                    sound_velocity_left_new = new_vson_value
                    self._function_to_vanish.eraseVariables()
            except ValueError as err:
                raise err
            # -----------------------------
            # Pression partie droite
            # -----------------------------
            density_current_value = self.density.current_right_value[self._enriched]
            density_new_value = self.density.new_right_value[self._enriched]
            pressure_current_value = self.pressure.current_right_value[self._enriched]
            pseudo_current_value = self.pseudo.current_right_value[self._enriched]
            energy_current_value = self.energy.current_right_value[self._enriched]
            energy_new_value = self.energy.new_right_value[self._enriched]
            shape = energy_new_value.shape
            nbr_cells_to_solve = shape[0]
            solution_value = np.zeros(shape, dtype=np.float64, order='C')
            new_pressure_value = np.zeros(shape, dtype=np.float64, order='C')
            new_vson_value = np.zeros(shape, dtype=np.float64, order='C')
            try:
                if self._external_library is not None:
                    pb_size = ctypes.c_int()
                    #
                    old_density = density_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    new_density = density_new_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    tmp = (pressure_current_value + 2. * pseudo_current_value)
                    pressure = tmp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    old_energy = energy_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    pb_size.value = nbr_cells_to_solve
                    solution = solution_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    new_pressure = new_pressure_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    new_vson = new_vson_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                    self._computePressureExternal(old_density, new_density, pressure, old_energy, pb_size,
                                                  solution, new_pressure, new_vson)
                    energy_right_new = solution[0:nbr_cells_to_solve]
                    pressure_right_new = new_pressure[0:nbr_cells_to_solve]
                    sound_velocity_right_new = new_vson[0:nbr_cells_to_solve]
                else:
                    my_variables = {'EquationOfState': self.proprietes.material.eos,
                                    'OldDensity': density_current_value,
                                    'NewDensity': density_new_value,
                                    'Pressure': pressure_current_value + 2. * pseudo_current_value,
                                    'OldEnergy': energy_current_value}
                    self._function_to_vanish.setVariables(my_variables)
                    solution = self._solver.computeSolution(energy_current_value)
                    new_pressure_value, _, new_vson_value = \
                        self.proprietes.material.eos.solveVolumeEnergy(1. / density_new_value, solution)
                    energy_right_new = solution
                    pressure_right_new = new_pressure_value
                    sound_velocity_right_new = new_vson_value
                    self._function_to_vanish.eraseVariables()
            except ValueError as err:
                raise err
            self.pressure.new_value[self._enriched] = \
                from_geometry_to_classic_field(pressure_left_new, pressure_right_new)
            self.pressure.new_enr_value[self._enriched] = \
                from_geometry_to_enrich_field(pressure_left_new, pressure_right_new)
            #
            self.energy.new_value[self._enriched] = \
                from_geometry_to_classic_field(energy_left_new, energy_right_new)
            self.energy.new_enr_value[self._enriched] = \
                from_geometry_to_enrich_field(energy_left_new, energy_right_new)
            #
            self.sound_velocity.new_value[self._enriched] = \
                from_geometry_to_classic_field(sound_velocity_left_new, sound_velocity_right_new)
            self.sound_velocity.new_enr_value[self._enriched] = \
                from_geometry_to_enrich_field(sound_velocity_left_new, sound_velocity_right_new)

    def computeNewSize(self, topologie, vecteur_coord_noeuds, vecteur_vitesse_noeuds,
                       vecteur_vitesse_enr_noeud, time_step = None):
        """
        Calcul des nouvelles longueurs de l'élément
        """
        # Calcul des tailles des éléments classiques
        super(OneDimensionEnrichedOneDimensionElement, self).computeNewSize(topologie, vecteur_coord_noeuds,
                                                                            mask=self._classical,
                                                                            time_step=time_step)
        if self._enriched.any():
            # Calcul des tailles des parties gauches des éléments enrichis
            connectivity = np.array(topologie._nodes_belonging_to_cell)[self._enriched]
            self.left_size.new_value[self._enriched] = self.left_size.current_value[self._enriched] + \
                                                       (0.5 * (vecteur_vitesse_noeuds[connectivity[:,1]] - vecteur_vitesse_enr_noeud[connectivity[:,0]]) -
                 0.5 * (vecteur_vitesse_noeuds[connectivity[:,0]] - vecteur_vitesse_enr_noeud[connectivity[:,0]]))\
                 * time_step
            self.right_size.new_value[self._enriched] = self.right_size.current_value[self._enriched] + \
                                                        (0.5 * (vecteur_vitesse_noeuds[connectivity[:,1]] - vecteur_vitesse_enr_noeud[connectivity[:,1]]) -
                 0.5 * (vecteur_vitesse_noeuds[connectivity[:,0]] - vecteur_vitesse_enr_noeud[connectivity[:,1]]))\
                 * time_step

    def computeNewDensity(self):
        """
        Calcul des nouvelles densités
        """
        #  Calcul des densités des éléments classiques
        super(OneDimensionEnrichedOneDimensionElement, self).computeNewDensity(mask=self._classical)
        #
        if self._enriched.any():
            densite_gauche_t_plus_dt = self.density.current_left_value[self._enriched] * \
                                       self.left_size.current_value[self._enriched] / self.left_size.new_value[
                                           self._enriched]
            densite_droite_t_plus_dt = self.density.current_right_value[self._enriched] * \
                                       self.right_size.current_value[self._enriched] / self.right_size.new_value[
                                           self._enriched]
            self.density.new_value[self._enriched] = \
                from_geometry_to_classic_field(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)
            self.density.new_enr_value[self._enriched] = \
                from_geometry_to_enrich_field(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)

    def computeNewPseudo(self, delta_t):
        """
        Calcul de la nouvelle pseudo
        """
        # calcul de la pseudo des éléments classique
        super(OneDimensionEnrichedOneDimensionElement, self).computeNewPseudo(delta_t, mask=self._classical)
        #
        if self._enriched.any():
            rho_t_gauche = self.density.current_left_value[self._enriched]
            rho_t_plus_dt_gauche = self.density.new_left_value[self._enriched]
            cson_t_gauche = self.sound_velocity.current_left_value[self._enriched]
            pseudo_gauche = \
                OneDimensionElement.computePseudo(delta_t, rho_t_gauche,
                                                  rho_t_plus_dt_gauche,
                                                  self.left_size.new_value[self._enriched],
                                                  cson_t_gauche,
                                                  self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)
    
            rho_t_droite = self.density.current_right_value[self._enriched]
            rho_t_plus_dt_droite = self.density.new_right_value[self._enriched]
            cson_t_droite = self.sound_velocity.current_right_value[self._enriched]
            pseudo_droite = \
                OneDimensionElement.computePseudo(delta_t, rho_t_droite,
                                                  rho_t_plus_dt_droite,
                                                  self.right_size.new_value[self._enriched],
                                                  cson_t_droite,
                                                  self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)
    
            self.pseudo.new_value[self._enriched] = \
                from_geometry_to_classic_field(pseudo_gauche, pseudo_droite)
            self.pseudo.new_enr_value[self._enriched] = \
                from_geometry_to_enrich_field(pseudo_gauche, pseudo_droite)

    def computeNewTimeStep(self):
        """
        Calcul du pas de temps
        """
        # calcul du pas de temps pour les éléments classiques
        super(OneDimensionEnrichedOneDimensionElement, self).computeNewTimeStep(mask=self._classical)
        #
        if self._enriched.any():
            cfl = self.proprietes.numeric.cfl
            rho_t_gauche = self.density.current_left_value[self._enriched]
            rho_t_plus_dt_gauche = self.density.new_left_value[self._enriched]
            cson_t_plus_dt_gauche = self.sound_velocity.new_left_value[self._enriched]
            pseudo_gauche = self.pseudo.current_left_value[self._enriched]
            dt_g = \
                OneDimensionElement.computeTimeStep(cfl, rho_t_gauche,
                                                    rho_t_plus_dt_gauche,
                                                    self.left_size.new_value[self._enriched],
                                                    cson_t_plus_dt_gauche,
                                                    pseudo_gauche)
    
            rho_t_droite = self.density.current_right_value[self._enriched]
            rho_t_plus_dt_droite = self.density.new_right_value[self._enriched]
            cson_t_plus_dt_droite = self.sound_velocity.new_right_value[self._enriched]
            pseudo_droite = self.pseudo.current_right_value[self._enriched]
            dt_d = \
                OneDimensionElement.computeTimeStep(cfl, rho_t_droite,
                                                    rho_t_plus_dt_droite,
                                                    self.right_size.new_value[self._enriched],
                                                    cson_t_plus_dt_droite,
                                                    pseudo_droite)

            self._dt[self._enriched] = dt_g + dt_d  # Bizarre --> A vérifier
