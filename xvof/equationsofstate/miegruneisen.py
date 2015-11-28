#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant une équation d'état de type Mie-Gruneisen
"""
import ctypes
import os
import numpy as np

from xvof.equationsofstate.equationofstatebase import EquationOfStateBase


# Deactivate pylint warnings due to NotImplementedError
# pylint: disable=R0921
class MieGruneisen(EquationOfStateBase):
    """
    Un objet décrivant l'équation d'état de type Mie_Gruneisen
    """
    def __init__(self, **kwargs):
        """
        self.czero    : paramètre czero
        self.coeff_s1       : paramètre S1
        self.coeff_s2       : paramètre S2
        self.coeff_s3       : paramètre S3
        self.rhozeo   : paramètre rhozero
        self.grunzero : paramètre grunzero
        self.coeff_b        : paramètre b
        self.ezero    : paramètre ezero
        """
        super(MieGruneisen, self).__init__()
        #
        self.__parameters = {'czero': 3980.0,
                             'S1': 1.58,
                             'S2': 0.,
                             'S3': 0.,
                             'rhozero': 8129.0,
                             'grunzero': 1.6,
                             'b': 0.5,
                             'ezero': 0.}
        for (prop, default) in self.__parameters.iteritems():
            self.__parameters[prop] = kwargs.get(prop, default)

#        _file = 'libMieGruneisen.so'
#        _path = os.path.join(*(os.path.split(__file__)[:-1] + (_file,)))
#        self._mod = ctypes.cdll.LoadLibrary(_path)
#        self._solveVE = self._mod.solveVolumeEnergy
#        self._solveVE.argtypes = ([ctypes.c_double, ] * 10 + [ctypes.POINTER(ctypes.c_double), ] * 3)

    @property
    def czero(self):
        """
        Vitesse du son standard
        """
        return self.__parameters['czero']

    @property
    def coeff_s1(self):
        """
        Paramètre S1
        """
        return self.__parameters['S1']

    @property
    def coeff_s2(self):
        """
        Paramètre S2
        """
        return self.__parameters['S2']

    @property
    def coeff_s3(self):
        """
        Paramètre S3
        """
        return self.__parameters['S3']

    @property
    def rhozero(self):
        """
        Masse volumique standard
        """
        return self.__parameters['rhozero']

    @property
    def grunzero(self):
        """
        Coefficient de Gruneisen standard
        """
        return self.__parameters['grunzero']

    @property
    def coeff_b(self):
        """
        Paramètre b
        """
        return self.__parameters['b']

    @property
    def ezero(self):
        """
        Energie interne standard
        """
        return self.__parameters['ezero']

    def __str__(self):
        message = "EquationOfState : {}".format(self.__class__)
        message += "\nParameters : "
        for key, val in sorted(self.__parameters.items(), key=lambda m: m[0]):
            message += "\n -- {:>20s} : {:>9.8g}".format(key, val)
        return message

    def __repr__(self):
        return self.__str__()

    def solveVolumeEnergy(self, specific_volume, internal_energy):
        """
        Fournit le triplet (pression | dérivée de la pression en
        fonction de l'énergie | vitesse du son) à partir du couple
        ( volume massique | énergie interne )
        """
        # -------------------------------------------------
        # Définition de variable locales pour eviter de
        # multiples recherche (gain de temps)
        # -------------------------------------------------
        rhozero = self.rhozero
        czero2 = self.czero ** 2
        # Dérivee du coefficient de gruneisen
        dgam = rhozero * (self.grunzero - self.coeff_b)
        #
        epsv = 1 - rhozero * specific_volume
        #
        epsv2 = epsv ** 2
        # Coefficient de gruneisen
        gampervol = (self.grunzero * (1 - epsv) + self.coeff_b * epsv) / specific_volume
        # -------------------------------------------------
        # Définition de variable locales redondantes
        # -------------------------------------------------
        redond_a = self.coeff_s1 + 2. * self.coeff_s2 * epsv + 3. * self.coeff_s3 * epsv2
        #
#         denom = np.zeros(specific_volume.shape, dtype=np.float64, order='C')
        phi = np.zeros(specific_volume.shape, dtype=np.float64, order='C')
        einth = np.zeros(specific_volume.shape, dtype=np.float64, order='C')
        dphi = np.zeros(specific_volume.shape, dtype=np.float64, order='C')
        deinth = np.zeros(specific_volume.shape, dtype=np.float64, order='C')
        dpdv = np.zeros(specific_volume.shape, dtype=np.float64, order='C')
        targets = epsv > 0
        #
#         if epsv > 0:
        # ============================================================
        # si epsv > 0, la pression depend de einth et phi.
        # einth : energie interne specifique sur l hugoniot
        # phi : pression sur l hugoniot
        # denom : racine du denominateur de phi
        # dphi : derivee de ph par rapport a
        # deinth : derivee de einth par rapport a v
        # dpdv : dp/dv
        # ============================================================
        denom = (1. - (self.coeff_s1 + self.coeff_s2 * epsv[targets] + self.coeff_s3 * epsv2[targets]) * epsv[targets])
        phi[targets] = rhozero * czero2 * epsv[targets] / denom ** 2
        einth[targets] = self.ezero + phi[targets] * epsv[targets] / (2. * rhozero)
        #
        dphi[targets] = phi[targets] * rhozero * (-1. / epsv[targets] - 2. * redond_a[targets] / denom)
        #
        deinth[targets] = phi[targets] * (-1. - epsv[targets] * redond_a[targets] / denom)
        #
        dpdv[targets] = dphi[targets] + (dgam - gampervol[targets]) * \
            (internal_energy[targets] - einth[targets]) / specific_volume[targets] - \
            gampervol[targets] * deinth[targets]
            #
#         elif epsv <= 0:
        targets = epsv <= 0
        # ============================================================
        # traitement en tension : epsv < 0
        # la pression depend d une fonction de v : phi
        # et
        # de e0 appelee einth
        # ============================================================
        phi[targets] = rhozero * czero2 * epsv[targets] / (1. - epsv[targets])
        # einth ---> e0
        einth[targets] = self.ezero
        #
        dphi[targets] = -czero2 / specific_volume[targets] ** 2
        #
        dpdv[targets] = dphi[targets] + (dgam - gampervol[targets]) * \
            (internal_energy[targets] - einth[targets]) / specific_volume[targets]
        # ****************************
        # Pression :
        # ****************************
        pressure = phi + (gampervol) * (internal_energy - einth)
        # ======================================
        # Carre de la vitesse du son :
        # ======================================
        vson2 = specific_volume ** 2 * (pressure * gampervol - dpdv)
        pb = vson2 < 0
        if pb.any():
            msg = "Carré de la vitesse du son < 0\n"
            msg += "specific_volume = {:15.9g}\n".format(specific_volume)
            msg += "pressure = {:15.9g}\n".format(pressure)
            msg += "dpsurde = {:15.9g}\n".format(gampervol)
            msg += "dpdv = {:15.9g}\n".format(dpdv)
            raise ValueError(msg)
        vson = vson2 ** 0.5
        #
        pb = np.logical_and((vson > 0.), (vson < 10000.))
        vson[~pb] = 0.
        #
        return pressure, gampervol, vson

#    def solveVolumeEnergy(self, specific_volume, internal_energy):
#        """
#        Fournit le triplet (pression | dérivée de la pression en
#        fonction de l'énergie | vitesse du son) à partir du couple
#        ( volume massique | énergie interne ) en passant par une lib C externe
#        """
#        pression = ctypes.c_double()
#        gamma_per_vol = ctypes.c_double()
#        vitson = ctypes.c_double()
#        self._solveVE(self.czero, self.coeff_s1, self.coeff_s2, self.coeff_s3,
#                      self.rhozero, self.grunzero, self.coeff_b, self.ezero,
#                      specific_volume, internal_energy, pression, gamma_per_vol,
#                      vitson)
#        return pression.value, gamma_per_vol.value, vitson.value

    def solveVolumePressure(self, specific_volume, pressure):
        raise NotImplementedError

    def solveVolumeTemperature(self, specific_volume, temperatue):
        raise NotImplementedError
