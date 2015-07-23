#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un élément en 1d
"""
# --------------------------------------------------------
#               IMPORTATIONS DIVERSES                    #
# --------------------------------------------------------
from xvof.element import Element
from xvof.solver.newtonraphson import NewtonRaphson
from xvof.solver.functionstosolve.vnrenergyevolutionforveformulation import VnrEnergyEvolutionForVolumeEnergyFormulation


# --------------------------------------------------------
#        DEFINITION DES CLASSES ET FONCTIONS             #
# --------------------------------------------------------
class Element1d(Element):
    """
    Une classe pour les éléments en 1D
    """
    nbr_noeuds = 2

    @classmethod
    def computePseudo(cls, delta_t, rho_old, rho_new, size_new,
                      cel_son, a_pseudo, b_pseudo):
        """
        Calcul de la pseudo
        """
        # pylint: disable=too-many-arguments
        # 8 arguments semblent nécessaires...
        vnt = 1. / rho_old
        vnplusun = 1. / rho_new
        vnplusundemi = 0.5 * (vnt + vnplusun)
        vpointnplusundemi = 1. / delta_t * (vnplusun - vnt)
        divu = vpointnplusundemi / vnplusundemi
        pseudo = 0.
        if divu < 0.:
            pseudo = 1. / vnplusundemi * \
                (a_pseudo * size_new ** 2 * vpointnplusundemi ** 2 / vnplusundemi ** 2 +
                 b_pseudo * size_new * cel_son *
                 abs(vpointnplusundemi) / vnplusundemi)
        return pseudo

    @classmethod
    def computeTimeStep(cls, cfl, rho_old, rho_new, taille_new, cson_new,
                        pseudo):
        """
        Calcul du pas de temps
        """
        # pylint: disable=too-many-arguments
        # 7 arguments pour cette méthode cela semble ok
        delta_t = 0.
        if (rho_new - rho_old) > 0.1:
            delta_t = cfl * taille_new / ((cson_new ** 2 + 2. * pseudo /
                                           (rho_new - rho_old)) ** 0.5)
        else:
            delta_t = cfl * taille_new / (cson_new)
        return delta_t

    def __init__(self, proprietes):
        Element.__init__(self, proprietes)
        self._function_to_vanish = VnrEnergyEvolutionForVolumeEnergyFormulation()
        self._solver = NewtonRaphson(self._function_to_vanish)

    # --------------------------------------------------------
    #            DEFINITION DES PROPRIETES                   #
    # --------------------------------------------------------
    @property
    def masse(self):
        """ Masse de l'élément """
        return self.taille_t * self.proprietes.geometric.section * \
            self.density.current_value

    # --------------------------------------------------------
    #            DEFINITION DES METHODES                     #
    # --------------------------------------------------------
    def computeNewPressure(self):
        """
        Calcul du triplet energie, pression, vitesse du son
        au pas de temps suivant
        Formulation v-e
        """
        try:
            my_variables = {'EquationOfState': self.proprietes.material.eos,
                            'OldDensity': self.density.current_value,
                            'NewDensity': self.density.new_value,
                            'Pressure': self.pressure.current_value + 2. * self.pseudo.current_value,
                            'OldEnergy': self.energy.current_value}
            self._function_to_vanish.setVariables(my_variables)
            self.energy.new_value = self._solver.computeSolution(self.energy.current_value)
            self.pressure.new_value, _, self.sound_velocity.new_value = \
                self.proprietes.material.eos.solveVolumeEnergy(1. / self.density.new_value, self.energy.new_value)
            self._function_to_vanish.eraseVariables()
        except ValueError as err:
            print "Element concerné : {}".format(self)
            raise err

    def computeSize(self, noeuds):
        """
        Calcul de la longueur de l'élément (à t)
        """
        self._size_t = abs(noeuds[0].coordt[0] -
                           noeuds[1].coordt[0])

    def computeNewSize(self, noeuds, time_step=None):
        """
        Calcul de la nouvelle longueur de l'élément (à t+dt)
        """
        self._size_t_plus_dt = abs(noeuds[0].coordtpdt[0] -
                                   noeuds[1].coordtpdt[0])

    def computeNewDensity(self):
        """
        Calcul de la densité à l'instant t+dt basé sur
        la conservation de la masse
        """
        self.density.new_value = \
            self.density.current_value * self.taille_t / self.taille_t_plus_dt

    def computeNewPseudo(self, delta_t):
        """
        Calcul de la nouvelle pseudo
        """
        self.pseudo.new_value = \
            Element1d.computePseudo(delta_t, self.density.current_value, self.density.new_value,
                                    self.taille_t_plus_dt, self.sound_velocity.current_value,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

    def computeNewTimeStep(self):
        """
        Calcul du pas de temps dans l'élément
        """
        cfl = self.proprietes.numeric.cfl
        self._dt = \
            Element1d.computeTimeStep(cfl, self.density.current_value, self.density.new_value,
                                      self.taille_t_plus_dt, self.sound_velocity.new_value,
                                      self.pseudo.current_value)

    def imposePressure(self, pression):
        """
        On impose la pression à t+dt (par exemple pour endommagement)
        """
        self.pressure.new_value = pression
