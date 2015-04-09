#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un élément en 1d
"""
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from xvof.element import Element
import numpy as np

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


class Element1d(Element):
    """
    Une classe pour les éléments en 1D
    """
    @classmethod
    def newton_raphson_for_ve(cls, eos, rho_old, rho_new, pression_old,
                              pseudo_old, nrj_old):
        """
        Algorithme de Newton-Raphson pour déterminer le couple
        energie/pression au pas de temps suivant
        Formulation v-e
        """
        delta_v = 1. / rho_new - 1. / rho_old
        pression_t = pression_old + 2. * pseudo_old
        # Variable du Newton
        nrj_i = nrj_old
        # Critère de convergence
        convergence = False
        # Nombre d'itérations
        nit = 0
        #

        def calcul_f_et_df(enerj):
            """
            Fonction à annuler et sa dérivée pour le schéma VNR
            Formulation v-e
            """
            (p_i, dpsurde, dummy) = eos.solve_ve(1. / rho_new, enerj)
            # Fonction à annuler
            func = enerj + p_i * delta_v / 2. + pression_t * delta_v / 2. -\
                nrj_old
            # Dérivée de la fonction à annuler
            dfunc = 1 + dpsurde * delta_v / 2.
            return (func, dfunc)
        #
        (func_i, dfunc_i_surde) = calcul_f_et_df(nrj_i)
        #
        while(not convergence and (nit < 100)):
            # Correction
            nrj_iplus1 = nrj_i - func_i / dfunc_i_surde
            nit += 1
            if(abs(func_i) < 1e-09):
                convergence = True
                res_nrj = nrj_i
                res_pression_t_plus_dt, dummy, res_cson = \
                eos.solve_ve(1. / rho_new, res_nrj)
                break
            # Incrémentation
            nrj_i = nrj_iplus1
            #
            (func_i, dfunc_i_surde) = calcul_f_et_df(nrj_i)
            if(abs(dfunc_i_surde) < 1.e-09):
                print "Sortie du NR par manque de pente :-)"
                convergence = True
                res_nrj = nrj_i
                res_pression_t_plus_dt, dummy, res_cson = \
                eos.solve_ve(1. / rho_new, res_nrj)
                break
        if(nit == 100):
            print "Erreur de convergence du NR"
            print "func_i=", func_i
            print "nit=", nit
            exit(255)
        return res_nrj, res_pression_t_plus_dt, res_cson

    def __init__(self, proprietes, indice, noeuds):
        Element.__init__(self, proprietes, indice, noeuds)
        self.noeuds = noeuds
        self._size_t = abs(self.noeuds[0].coordtpdt[0] -
            self.noeuds[1].coordtpdt[0])

    #------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    #------------------------------------------------------------
    @property
    def noeuds(self):
        """
        Getter de la liste des noeuds de l'élément
        """
        return self._noeuds

    @noeuds.setter
    def noeuds(self, list_noeuds):
        """
        Setter de la liste des noeuds de l'élément
        """
        if (len(list_noeuds) != 2):
            raise SystemExit("En 1D, un élément possède 2 noeuds!")
        self._noeuds[:] = list_noeuds[:]
        self._noeuds = \
        sorted(self._noeuds, key=lambda m: m.coordt[0])

    #------------------------------------------------------------
    # DEFINITIONS DES METHODES
    #------------------------------------------------------------
    def calculer_nouvo_pression(self):
        """
        Calcul du triplet energie, pression, vitesse du son
        au pas de temps suivant
        Formulation v-e

        TEST UNITAIRE
        >>> from xvof.miscellaneous import *
        >>> from xvof.equationsofstate import MieGruneisen
        >>> ee = MieGruneisen()
        >>> num_props = numerical_props(0.2, 1.0, 0.35)
        >>> mat_props = material_props(1.0e+05, 0.0, 8129., ee)
        >>> geom_props = geometrical_props(1.0e-06)
        >>> props = properties(num_props, mat_props, geom_props)
        >>> my_elem = Element1d(props, 123, [])
        >>> my_elem._rho_t_plus_dt = 9000.0
        >>> my_elem.calculer_nouvo_pression()
        >>> print my_elem.nrj_t_plus_dt/1e+05
        1.03378424404
        >>> print my_elem.pression_t_plus_dt/1.0e+09
        17.3667631638
        >>> print my_elem.cson_t_plus_dt/1.0e+03
        4.89803134405
        """
        self._nrj_t_plus_dt, self._pression_t_plus_dt, self._cson_t_plus_dt = \
            Element1d.newton_raphson_for_ve(self.proprietes.material.eos,
                self.rho_t, self.rho_t_plus_dt, self.pression_t, self.pseudo,
                self.nrj_t)

    def calculer_nouvo_taille(self):
        """
        Calcul de la nouvelle longueur de l'élément

        TEST UNITAIRE
        >>> from xvof.miscellaneous import *
        >>> from xvof.equationsofstate import MieGruneisen
        >>> ee = MieGruneisen()
        >>> num_props = numerical_props(0.2, 1.0, 0.35)
        >>> mat_props = material_props(1.0e+05, 0.0, 8129., ee)
        >>> geom_props = geometrical_props(1.0e-06)
        >>> props = properties(num_props, mat_props, geom_props)
        >>> my_elem = Element1d(props, 123, 2.5e-03)
        >>> class noeuds:
        ...     pass
        >>> noe1 = noeuds()
        >>> noe1.coordtpdt = np.array([5.0e-03])
        >>> noe2 = noeuds()
        >>> noe2.coordtpdt = np.array([6.0e-03])
        >>> my_elem.noeuds = [noe2, noe1]
        >>> my_elem.calculer_nouvo_taille()
        >>> print my_elem.taille_t_plus_dt
        0.001
        """
        self._size_t_plus_dt = abs(self.noeuds[0].coordtpdt[0] -
                                    self.noeuds[1].coordtpdt[0])

    def calculer_nouvo_densite(self):
        """
        Calcul de la densité à l'instant t+dt basé sur
        la conservation de la masse

        TEST UNITAIRE
        >>> import numpy as np
        >>> from xvof.node import Node1d
        >>> from xvof.miscellaneous import *
        >>> from xvof.equationsofstate import MieGruneisen
        >>> ee = MieGruneisen()
        >>> num_props = numerical_props(0.2, 1.0, 0.35)
        >>> mat_props = material_props(1.0e+05, 0.0, 8129., ee)
        >>> geom_props = geometrical_props(1.0e-06)
        >>> props = properties(num_props, mat_props, geom_props)
        >>> noda = Node1d(1, poz_init=np.array([-1.0e-03]))
        >>> nodb = Node1d(2, poz_init=np.array([1.5e-03]))
        >>> my_elem = Element1d(props, 123, [nodb, noda])
        >>> my_elem._size_t_plus_dt = 1.25e-03
        >>> my_elem.calculer_nouvo_densite()
        >>> print my_elem.rho_t_plus_dt
        16258.0
        """
        self._rho_t_plus_dt = \
            self.rho_t * self.taille_t / self.taille_t_plus_dt

    def calculer_nouvo_pseudo(self, delta_t):
        """
        Calcul de la nouvelle pseudo

        TEST UNITAIRE
        >>> import numpy as np
        >>> from xvof.node import Node1d
        >>> from xvof.miscellaneous import *
        >>> from xvof.equationsofstate import MieGruneisen
        >>> ee = MieGruneisen()
        >>> num_props = numerical_props(0.2, 1.0, 0.35)
        >>> mat_props = material_props(1.0e+05, 0.0, 8129., ee)
        >>> geom_props = geometrical_props(1.0e-06)
        >>> props = properties(num_props, mat_props, geom_props)
        >>> noe1 = Node1d(1, poz_init=np.array([4.0e-03]))
        >>> noe2 = Node1d(2, poz_init=np.array([7.0e-03]))
        >>> noe1.coordtpdt = np.array([5.0e-03])
        >>> noe2.coordtpdt = np.array([6.0e-03])
        >>> my_elem = Element1d(props, 123, [noe2, noe1])
        >>> my_elem.calculer_nouvo_taille()
        >>> my_elem.calculer_nouvo_densite()
        >>> my_elem.calculer_nouvo_pseudo(1.0e-6)
        >>> print my_elem.pseudo
        1706379008.75
        """
        vnt = 1. / self.rho_t
        vnplusun = 1. / self.rho_t_plus_dt
        vnplusundemi = 0.5 * (vnt + vnplusun)
        vpointnplusundemi = 1. / delta_t * (vnplusun - vnt)
        divu = vpointnplusundemi / vnplusundemi
        pseuda = self.proprietes.numeric.a_pseudo
        pseudb = self.proprietes.numeric.b_pseudo
        if(divu < 0.):
            self._pseudo_plus_un_demi = 1. / vnplusundemi * \
            (
            pseuda * self.taille_t_plus_dt ** 2 * vpointnplusundemi ** 2 /
            vnplusundemi ** 2 +
            pseudb * self.taille_t_plus_dt * self.cson_t *
            abs(vpointnplusundemi) / vnplusundemi
            )
        else:
            self._pseudo_plus_undemi = 0.

    def calculer_nouvo_dt(self):
        """
        Calcul du pas de temps dans l'élément

        TEST UNITAIRE
        >>> import numpy as np
        >>> from xvof.node import Node1d
        >>> from xvof.miscellaneous import *
        >>> from xvof.equationsofstate import MieGruneisen
        >>> ee = MieGruneisen()
        >>> num_props = numerical_props(0.2, 1.0, 0.35)
        >>> mat_props = material_props(1.0e+05, 0.0, 8129., ee)
        >>> geom_props = geometrical_props(1.0e-06)
        >>> props = properties(num_props, mat_props, geom_props)
        >>> noe1 = Node1d(1, poz_init=np.array([4.0e-03]))
        >>> noe2 = Node1d(2, poz_init=np.array([7.0e-03]))
        >>> noe1._x_t_plus_dt = np.array([5.0e-03])
        >>> noe2._x_t_plus_dt = np.array([7.0e-03])
        >>> my_elem = Element1d(props, 123, [noe1, noe2])
        >>> my_elem.calculer_nouvo_taille()
        >>> my_elem.calculer_nouvo_densite()
        >>> my_elem.calculer_nouvo_pseudo(1.0e-6)
        >>> my_elem.calculer_nouvo_pression()
        >>> my_elem.calculer_nouvo_dt()
        >>> print my_elem.delta_t
        2.63819095477e-07
        """
        cfl = self.proprietes.numeric.cfl
        if((self.rho_t_plus_dt - self.rho_t) > 0.1):
            self._dt = cfl * self.taille_t_plus_dt / \
            ((self.cson_t_plus_dt ** 2 + 2. * self.pseudo /
            (self.rho_t_plus_dt - self.rho_t)) ** 0.5)
        else:
            self._dt = cfl * self.taille_t_plus_dt / (self.cson_t_plus_dt)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#######          PROGRAMME PRINCIPAL        ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if __name__ == "__main__":
    import doctest
    TESTRES = doctest.testmod(verbose=0)
    if(TESTRES[0] == 0):
        print "TESTS UNITAIRES : OK"
        #
        # Lancemet du profiling
        #
        from xvof.utilities import timeit_file
        from os import system
        from xvof.miscellaneous import numerical_props, material_props
        from xvof.miscellaneous import geometrical_props, properties
        from xvof.equationsofstate import MieGruneisen
        #

        @timeit_file('calculer_nouvo_pression.log')
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
        MY_ELEM = Element1d(PROPS, 123, 2.5e-03)
        MY_ELEM._rho_t_plus_dt = 9000.0
        NBR_ITER = 100000
        print "Lancement profiling sur {:6d} itérations".format(NBR_ITER)
        profil_calculer_nouvo_pression(MY_ELEM, NBR_ITER)
        system('cat calculer_nouvo_pression.log')