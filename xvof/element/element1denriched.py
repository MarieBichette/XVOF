#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un élément enrichi en 1d
"""
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ########### IMPORTATIONS DIVERSES  ####################
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import numpy as np
from xvof.element import Element1d
from xvof.solver.newtonraphson import NewtonRaphson
from xvof.solver.functionstosolve.vnrenergyevolutionforveformulation import VnrEnergyEvolutionForVolumeEnergyFormulation


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ###### DEFINITION DES CLASSES & FONCTIONS  ###############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class Element1dEnriched(Element1d):
    """
    Une classe pour les éléments enrichis dans le cas 1d
    """
    @classmethod
    def fromGeometryToEnrichField(cls, champ_gauche, champ_droite):
        """
        Renvoi le champ enrichi à partir des champs gauche et droite
        """
        return (champ_droite - champ_gauche) * 0.5

    @classmethod
    def fromGeometryToClassicField(cls, champ_gauche, champ_droite):
        """
        Renvoi le champ classique à partir des champs gauche et droite
        """
        return (champ_droite + champ_gauche) * 0.5

    @classmethod
    def fromEnrichToLeftPartField(cls, champ_classic, champ_enrich):
        """
        Renvoi le champ à gauche d'après les champs classsique et enrichis
        """
        return champ_classic - champ_enrich

    @classmethod
    def fromEnrichToRightPartField(cls, champ_classic, champ_enrich):
        """
        Renvoi le champ à droite d'après les champs classsique et enrichis
        """
        return champ_classic + champ_enrich

    def getLeftField(self, key_field):
        fields = {'Pressure': (self._pression_t, self._pression_t_enrichi),
                  'NewPressure': (self._pression_t_plus_dt, self._pression_t_plus_dt_enrichi),
                  'Density': (self._rho_t, self._rho_t_enrichi),
                  'NewDensity': (self._rho_t_plus_dt, self._rho_t_plus_dt_enrichi),
                  'Energy': (self._nrj_t, self._nrj_t_enrichi),
                  'Pseudo': (self._pseudo_plus_un_demi, self._pseudo_plus_un_demi_enrichi),
                  'SoundVelocity': (self._cson_t, self._cson_t_enrichi),
                  'NewSoundVelocity': (self._cson_t_plus_dt, self._cson_t_plus_dt_enrichi)}
        return Element1dEnriched.fromEnrichToLeftPartField(fields[key_field][0], fields[key_field][1])

    def getRightField(self, key_field):
        fields = {'Pressure': (self._pression_t, self._pression_t_enrichi),
                  'NewPressure': (self._pression_t_plus_dt, self._pression_t_plus_dt_enrichi),
                  'Density': (self._rho_t, self._rho_t_enrichi),
                  'NewDensity': (self._rho_t_plus_dt, self._rho_t_plus_dt_enrichi),
                  'Energy': (self._nrj_t, self._nrj_t_enrichi),
                  'Pseudo': (self._pseudo_plus_un_demi, self._pseudo_plus_un_demi_enrichi),
                  'SoundVelocity': (self._cson_t, self._cson_t_enrichi),
                  'NewSoundVelocity': (self._cson_t_plus_dt, self._cson_t_plus_dt_enrichi)}
        return Element1dEnriched.fromEnrichToRightPartField(fields[key_field][0], fields[key_field][1])

    def __init__(self, element_origin, pos_discontin):
        Element1d.__init__(self, element_origin.proprietes)
        self._function_to_vanish = VnrEnergyEvolutionForVolumeEnergyFormulation()
        self._solver = NewtonRaphson(self._function_to_vanish, self.nrj_t)
        #
        if(pos_discontin < 0.) or (pos_discontin > 1.):
            message = "La position de la discontinuité dans"
            message += " l'élément enrichi doit être comprise entre 0 et 1!"
            raise SystemExit(message)
        #
        self._index = element_origin.index
        self._size_t = element_origin.taille_t
        self._size_t_plus_dt = element_origin.taille_t_plus_dt
        self._pression_t = element_origin.pression_t
        self._pression_t_plus_dt = element_origin.pression_t_plus_dt
        self._rho_t = element_origin.rho_t
        self._rho_t_plus_dt = element_origin.rho_t_plus_dt
        self._nrj_t = element_origin.nrj_t
        self._nrj_t_plus_dt = element_origin.nrj_t_plus_dt
        self._pseudo_plus_un_demi = element_origin.pseudo
        self._cson_t = element_origin.cson_t
        self._cson_t_plus_dt = element_origin.cson_t_plus_dt
        #
        self._pression_t_enrichi = 0.
        self._pression_t_plus_dt_enrichi = 0.
        self._rho_t_enrichi = 0.
        self._rho_t_plus_dt_enrichi = 0.
        self._nrj_t_enrichi = 0.
        self._nrj_t_plus_dt_enrichi = 0.
        self._pseudo_plus_un_demi_enrichi = 0.
        self._cson_t_enrichi = 0.
        self._cson_t_plus_dt_enrichi = 0.
        #
        self._taille_gauche_t = element_origin.taille_t * pos_discontin
        self._taille_gauche_t_plus_dt = element_origin.taille_t_plus_dt * pos_discontin
        self._taille_droite_t = element_origin.taille_t * (1. - pos_discontin)
        self._taille_droite_t_plus_dt = element_origin.taille_t_plus_dt * (1. - pos_discontin)
        #

    def getLeftPartCoordinates(self, noeuds):
        """
        Position du centre de l'élément au temps t
        """
        vec_coord = np.zeros(noeuds[0].dimension)
        vec_coord = noeuds[0].coordt[:] + self._taille_gauche_t / 2.0
        return vec_coord

    def getRightPartCoordinates(self, noeuds):
        """
        Position du centre de l'élément au temps t
        """
        vec_coord = np.zeros(noeuds[0].dimension)
        vec_coord = noeuds[1].coordt[:] - self._taille_droite_t / 2.0
        return vec_coord

    # --------------------------------------------------------
    #        DEFINITION DES METHODES                         #
    # --------------------------------------------------------
    def __str__(self):
        message = "ELEMENT ENRICHI {:4d} ".format(self._index)
        return message

    def printInfos(self):
        """
        Affichage des informations concernant l'élément
        """
        Element1d.printInfos(self)
        message = "==> masse volumique à gauche à t = {}\n".\
            format(Element1dEnriched.fromEnrichToLeftPartField(self._rho_t, self._rho_t_enrichi))
        message += "==> masse volumique à droite à t = {}\n".\
            format(Element1dEnriched.fromEnrichToRightPartField(self._rho_t, self._rho_t_enrichi))
        message += "==> masse volumique à gauche à t+dt = {}\n".\
            format(Element1dEnriched.fromEnrichToLeftPartField(self._rho_t_plus_dt, self._rho_t_plus_dt_enrichi))
        message += "==> masse volumique à droite à t+dt = {}\n".\
            format(Element1dEnriched.fromEnrichToRightPartField(self._rho_t_plus_dt, self._rho_t_plus_dt_enrichi))
        message += "==> taille à gauche à t = {}\n".\
            format(self._taille_gauche_t)
        message += "==> taille à droite à t = {}\n".\
            format(self._taille_droite_t)
        message += "==> taille à gauche à t+dt = {}\n".\
            format(self._taille_gauche_t_plus_dt)
        message += "==> taille à droite à t+dt = {}\n".\
            format(self._taille_droite_t_plus_dt)
        message += "==> pression à gauche à t = {}\n".\
            format(Element1dEnriched.fromEnrichToLeftPartField(self._pression_t, self._pression_t_enrichi))
        message += "==> pression à droite à t = {}\n".\
            format(Element1dEnriched.fromEnrichToRightPartField(self._pression_t, self._pression_t_enrichi))
        message += "==> vitesse du son à gauche à t = {}\n".\
            format(Element1dEnriched.fromEnrichToLeftPartField(self._cson_t, self._cson_t_enrichi))
        message += "==> vitesse du son à droite à t = {}\n".\
            format(Element1dEnriched.fromEnrichToRightPartField(self._cson_t, self._cson_t_enrichi))
        message += "==> vitesse du son à gauche à t+dt = {}\n".\
            format(Element1dEnriched.fromEnrichToLeftPartField(self._cson_t_plus_dt, self._cson_t_plus_dt_enrichi))
        message += "==> vitesse du son à droite à t+dt = {}\n".\
            format(Element1dEnriched.fromEnrichToRightPartField(self._cson_t_plus_dt, self._cson_t_plus_dt_enrichi))
        message += "==> énergie à gauche à t = {}\n".\
            format(Element1dEnriched.fromEnrichToLeftPartField(self._nrj_t, self._nrj_t_enrichi))
        message += "==> énergie à droite à t = {}\n".\
            format(Element1dEnriched.fromEnrichToRightPartField(self._nrj_t, self._nrj_t_enrichi))
        message += "==> énergie à gauche à t+dt = {}\n".\
            format(Element1dEnriched.fromEnrichToLeftPartField(self._nrj_t_plus_dt, self._nrj_t_plus_dt_enrichi))
        message += "==> énergie à droite à t+dt = {}\n".\
            format(Element1dEnriched.fromEnrichToRightPartField(self._nrj_t_plus_dt, self._nrj_t_plus_dt_enrichi))
        message += "==> pseudo à gauche = {}\n".\
            format(Element1dEnriched.fromEnrichToLeftPartField(self._pseudo_plus_un_demi,
                                                               self._pseudo_plus_un_demi_enrichi))
        message += "==> pseudo à droite = {}\n".\
            format(Element1dEnriched.fromEnrichToRightPartField(self._pseudo_plus_un_demi,
                                                                self._pseudo_plus_un_demi_enrichi))
        print message

    def computeNewPressure(self):
        """
        Calcul du triplet energie, pression, vitesse du son
        au pas de temps suivant
        Formulation v-e
        """
        # Traitement partie gauche
        my_variables = {'EquationOfState': self.proprietes.material.eos,
                        'OldDensity': self.getLeftField('Density'),
                        'NewDensity': self.getLeftField('NewDensity'),
                        'Pressure': self.getLeftField('Pressure') + 2. * self.getLeftField('Pseudo'),
                        'OldEnergy': self.getLeftField('Energy')}
        self._function_to_vanish.setVariables(my_variables)
        nrj_t_plus_dt_g = self._solver.computeSolution()
        pression_t_plus_dt_g, _, cson_t_plus_dt_g = \
            self.proprietes.material.eos.solveVolumeEnergy(1. / self.getLeftField('NewDensity'), nrj_t_plus_dt_g)
        self._function_to_vanish.eraseVariables()
        # Traitement partie droite
        my_variables = {'EquationOfState': self.proprietes.material.eos,
                        'OldDensity': self.getRightField('Density'),
                        'NewDensity': self.getRightField('NewDensity'),
                        'Pressure': self.getRightField('Pressure') + 2. * self.getRightField('Pseudo'),
                        'OldEnergy': self.getRightField('Energy')}
        self._function_to_vanish.setVariables(my_variables)
        nrj_t_plus_dt_d = self._solver.computeSolution()
        pression_t_plus_dt_d, _, cson_t_plus_dt_d = \
            self.proprietes.material.eos.solveVolumeEnergy(1. / self.getRightField('NewDensity'), nrj_t_plus_dt_d)
        self._function_to_vanish.eraseVariables()
        #
        self._pression_t_plus_dt = \
            Element1dEnriched.fromGeometryToClassicField(pression_t_plus_dt_g, pression_t_plus_dt_d)
        self._pression_t_plus_dt_enrichi = \
            Element1dEnriched.fromGeometryToEnrichField(pression_t_plus_dt_g, pression_t_plus_dt_d)
        #
        self._nrj_t_plus_dt = \
            Element1dEnriched.fromGeometryToClassicField(nrj_t_plus_dt_g, nrj_t_plus_dt_d)
        self._nrj_t_plus_dt_enrichi = \
            Element1dEnriched.fromGeometryToEnrichField(nrj_t_plus_dt_g, nrj_t_plus_dt_d)
        #
        self._cson_t_plus_dt = \
            Element1dEnriched.fromGeometryToClassicField(cson_t_plus_dt_g, cson_t_plus_dt_d)
        self._cson_t_plus_dt_enrichi = \
            Element1dEnriched.fromGeometryToEnrichField(cson_t_plus_dt_g, cson_t_plus_dt_d)

    def computeNewSize(self, noeuds, time_step=None):
        """
        Calcul des nouvelles longueurs de l'élément
        """
        # Les noeuds sont classés par coord croissante
        nod_g = noeuds[0]
        nod_d = noeuds[1]
        self._taille_gauche_t_plus_dt = self._taille_gauche_t + \
            (0.5 * (nod_d.upundemi_classique - nod_g.upundemi_enrichi) -
             0.5 * (nod_g.upundemi_classique - nod_g.upundemi_enrichi)) \
            * time_step
        self._taille_droite_t_plus_dt = self._taille_droite_t + \
            (0.5 * (nod_d.upundemi_classique + nod_d.upundemi_enrichi) -
             0.5 * (nod_g.upundemi_classique + nod_d.upundemi_enrichi)) \
            * time_step

    def computeNewDensity(self):
        """
        Calcul des nouvelles densités
        """
        densite_gauche_t_plus_dt = self.getLeftField('Density') * self._taille_gauche_t / self._taille_gauche_t_plus_dt
        densite_droite_t_plus_dt = self.getRightField('Density') * self._taille_droite_t / self._taille_droite_t_plus_dt
        self._rho_t_plus_dt = \
            Element1dEnriched.fromGeometryToClassicField(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)
        self._rho_t_plus_dt_enrichi = \
            Element1dEnriched.fromGeometryToEnrichField(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)

    def computeNewPseudo(self, delta_t):
        """
        Calcul de la nouvelle pseudo
        """
        rho_t_gauche = self.getLeftField('Density')
        rho_t_plus_dt_gauche = self.getLeftField('NewDensity')
        cson_t_gauche = self.getLeftField('SoundVelocity')
        pseudo_gauche = \
            Element1d.computePseudo(delta_t, rho_t_gauche,
                                    rho_t_plus_dt_gauche,
                                    self._taille_gauche_t_plus_dt,
                                    cson_t_gauche,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

        rho_t_droite = self.getRightField('Density')
        rho_t_plus_dt_droite = self.getRightField('NewDensity')
        cson_t_droite = self.getRightField('SoundVelocity')
        pseudo_droite = \
            Element1d.computePseudo(delta_t, rho_t_droite,
                                    rho_t_plus_dt_droite,
                                    self._taille_droite_t_plus_dt,
                                    cson_t_droite,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

        self._pseudo_plus_un_demi = \
            Element1dEnriched.fromGeometryToClassicField(pseudo_gauche, pseudo_droite)
        self._pseudo_plus_un_demi_enrichi = \
            Element1dEnriched.fromGeometryToEnrichField(pseudo_gauche, pseudo_droite)

    def computeNewTimeStep(self):
        """
        Calcul du pas de temps
        """
        cfl = self.proprietes.numeric.cfl
        rho_t_gauche = self.getLeftField('Density')
        rho_t_plus_dt_gauche = self.getLeftField('NewDensity')
        cson_t_plus_dt_gauche = self.getLeftField('NewSoundVelocity')
        pseudo_gauche = self.getLeftField('Pseudo')
        dt_g = \
            Element1d.computeTimeStep(cfl, rho_t_gauche,
                                      rho_t_plus_dt_gauche,
                                      self._taille_gauche_t_plus_dt,
                                      cson_t_plus_dt_gauche,
                                      pseudo_gauche)

        rho_t_droite = self.getRightField('Density')
        rho_t_plus_dt_droite = self.getRightField('NewDensity')
        cson_t_plus_dt_droite = self.getRightField('NewSoundVelocity')
        pseudo_droite = self.getRightField('Pseudo')
        dt_d = \
            Element1d.computeTimeStep(cfl, rho_t_droite,
                                      rho_t_plus_dt_droite,
                                      self._taille_droite_t_plus_dt,
                                      cson_t_plus_dt_droite,
                                      pseudo_droite)

        self._dt = dt_g + dt_d  # Bizarre --> A vérifier

    def incrementVariables(self):
        """
        Incrémentation des variables
        """
        Element1d.incrementVariables(self)
        self._pression_t_enrichi = self._pression_t_plus_dt_enrichi
        self._rho_t_enrichi = self._rho_t_plus_dt_enrichi
        self._cson_t_enrichi = self._cson_t_plus_dt_enrichi
        self._nrj_t_enrichi = self._nrj_t_plus_dt_enrichi
        self._taille_gauche_t = self._taille_gauche_t_plus_dt
        self._taille_droite_t = self._taille_droite_t_plus_dt
