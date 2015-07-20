#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un élément en 1d
"""
# --------------------------------------------------------
#               IMPORTATIONS DIVERSES                    #
# --------------------------------------------------------
import numpy as np
from xvof.element import Element


# --------------------------------------------------------
#        DEFINITION DES CLASSES ET FONCTIONS             #
# --------------------------------------------------------
class Element1d(Element):
    """
    Une classe pour les éléments en 1D
    """
    nbr_noeuds = 2

    @classmethod
    def computeFunctionAndDerivative(cls, eos, enerj, rho_new, rho_old, pression_t, enerj_old):
        """
        Fonction à annuler et sa dérivée pour le schéma VNR
        Formulation v-e
        """
        (p_i, dpsurde, dummy) = eos.solveVolumeEnergy(1. / rho_new, enerj)
        # Fonction à annuler
        delta_v = 1. / rho_new - 1. / rho_old
        func = enerj + p_i * delta_v / 2. + pression_t * delta_v / 2. - \
            enerj_old
        # Dérivée de la fonction à annuler
        dfunc = 1 + dpsurde * delta_v / 2.
        return (func, dfunc)

    @classmethod
    def executeNewtonRaphsonForVolumeEnergyFormulation(cls, function, eos, old_state, rho_new):
        """
        Algorithme de Newton-Raphson pour déterminer le couple
        energie/pression au pas de temps suivant
        Formulation v-e
        """
        rho_old = old_state['Density']
        pression_old = old_state['Pressure']
        pseudo_old = old_state['Pseudo']
        nrj_old = old_state['Energy']
        iter_max = 20
        pression_t = pression_old + 2. * pseudo_old
        # Variable du Newton
        nrj_i = nrj_old
        # Critère de convergence
        convergence = False
        # Nombre d'itérations
        nit = 0
        #
        while not convergence and (nit < iter_max):
            (func_i, dfunc_i_surde) = \
                function(eos, nrj_i, rho_new, rho_old, pression_t, nrj_old)
            # Correction
            nrj_iplus1 = nrj_i - func_i / dfunc_i_surde
            nit += 1
            if abs(func_i) < 1e-09:
                convergence = True
                res_nrj = nrj_i
                res_pression_t_plus_dt, dummy, res_cson = \
                    eos.solveVolumeEnergy(1. / rho_new, res_nrj)
                break
            if abs(dfunc_i_surde) < 1.e-09:
                print "Sortie du NR par manque de pente :-)"
                convergence = True
                res_nrj = nrj_i
                res_pression_t_plus_dt, dummy, res_cson = \
                    eos.solveVolumeEnergy(1. / rho_new, res_nrj)
                break
            # Incrémentation
            nrj_i = nrj_iplus1
        if nit == iter_max:
            print "Erreur de convergence du NR"
            print "func_i=", func_i
            print "nit=", nit
            raise ValueError()
        return res_nrj, res_pression_t_plus_dt, res_cson

    @classmethod
    def computePseudo(cls, delta_t, rho_old, rho_new, size_new,
                      cel_son, a_pseudo, b_pseudo):
        """
        Calcul de la pseudo
        """
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

    # --------------------------------------------------------
    #            DEFINITION DES PROPRIETES                   #
    # --------------------------------------------------------
    @property
    def masse(self):
        """ Masse de l'élément """
        return self.taille_t * self.proprietes.geometric.section * \
            self.rho_t

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
            eos = self.proprietes.material.eos
            function_to_vanish = Element1d.computeFunctionAndDerivative
            old_state = {'Pressure': self.pression_t,
                         'Density': self.rho_t,
                         'Pseudo': self.pseudo,
                         'Energy': self.nrj_t}
            self._nrj_t_plus_dt, self._pression_t_plus_dt, self._cson_t_plus_dt = \
                Element1d.executeNewtonRaphsonForVolumeEnergyFormulation(function_to_vanish, eos,
                                                                         old_state, self.rho_t_plus_dt)
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
        self._rho_t_plus_dt = \
            self.rho_t * self.taille_t / self.taille_t_plus_dt

    def computeNewPseudo(self, delta_t):
        """
        Calcul de la nouvelle pseudo
        """
        self._pseudo_plus_un_demi = \
            Element1d.computePseudo(delta_t, self.rho_t, self.rho_t_plus_dt,
                                    self.taille_t_plus_dt, self.cson_t,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

    def computeNewTimeStep(self):
        """
        Calcul du pas de temps dans l'élément
        """
        cfl = self.proprietes.numeric.cfl
        self._dt = \
            Element1d.computeTimeStep(cfl, self.rho_t, self.rho_t_plus_dt,
                                      self.taille_t_plus_dt, self.cson_t_plus_dt,
                                      self.pseudo)

    def imposePressure(self, pression):
        """
        On impose la pression à t+dt (par exemple pour endommagement)
        """
        self._pression_t_plus_dt = pression

# --------------------------------------------------------
#            PROGRAMME PRINCIPAL                         #
# --------------------------------------------------------
if __name__ == "__main__":
    #
    # Lancemet du profiling
    #
    from xvof.utilities import timeit_file
    from os import system
    from xvof.miscellaneous import numerical_props, material_props
    from xvof.miscellaneous import geometrical_props, properties
    from xvof.equationsofstate import MieGruneisen
    from xvof.node import Node1d
    #

    @timeit_file('/tmp/calculer_nouvo_pression.log')
    def profil_calculer_nouvo_pression(element, it_nb=100):
        """
        Fait it_nb appel(s) à calculer_nouvo_pression à
        des fins de profiling
        """
        i = 0
        while i < it_nb:
            element.calculer_nouvo_pression()
            i += 1
    #
    EE = MieGruneisen()
    NUM_PROPS = numerical_props(0.2, 1.0, 0.35)
    MAT_PROPS = material_props(1.0e+05, 0.0, 8129., EE)
    GEOM_PROPS = geometrical_props(1.0e-06)
    PROPS = properties(NUM_PROPS, MAT_PROPS, GEOM_PROPS)
    NODA = Node1d(1, np.array([0.5e-03]))
    NODB = Node1d(2, np.array([3.0e-03]))
    MY_ELEM = Element1d(PROPS)
    MY_ELEM._rho_t_plus_dt = 9000.0
    NBR_ITER = 100000
    print "Lancement profiling sur {:6d} itérations".format(NBR_ITER)
    profil_calculer_nouvo_pression(MY_ELEM, NBR_ITER)
    system('cat /tmp/calculer_nouvo_pression.log')
