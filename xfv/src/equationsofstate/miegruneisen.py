# -*- coding: iso-8859-1 -*-
"""
Implementing MieGruneisen class

>>> import numpy as np
>>> from xfv.src.equationsofstate.miegruneisen import MieGruneisen
>>> my_eos = MieGruneisen(3980, 1.58, 0, 0, 8129, 1.6, 0.5, 0)
>>> print my_eos # doctest:+NORMALIZE_WHITESPACE
EquationOfState : MieGruneisen
Parameters :
 --                czero :      3980
 --                   S1 :      1.58
 --                   S2 :         0
 --                   S3 :         0
 --              rhozero :      8129
 --             grunzero :       1.6
 --                    b :       0.5
 --                ezero :         0
 >>> my_eos # doctest:+NORMALIZE_WHITESPACE
 MieGruneisen(czero=3980.000000, S1=1.580000, S2=0.000000, S3=0.000000, rhozero=8129.000000, grunzero=1.600000,
 b=0.500000, ezero=0.000000)
 >>> density = np.array([9000., 8500., 9500.], dtype=np.float64, order='C')
 >>> specific_volume = 1. / density
 >>> internal_energy = np.array([1.0e+04, 1.0e+03, 1.0e+05], dtype=np.float64, order='C')
 >>> the_shape = internal_energy.shape
 >>> pressure = np.zeros(the_shape, dtype=np.float64, order='C')
 >>> sound_speed = np.zeros(the_shape, dtype=np.float64, order='C')
 >>> dpde = np.zeros(the_shape, dtype=np.float64, order='C')
 >>> my_eos.solveVolumeEnergy(specific_volume, internal_energy, pressure, sound_speed, dpde)
 >>> print pressure
 [  1.61115797e+10   6.26727977e+09   2.87613980e+10]
 >>> print sound_speed
 [ 4871.9323597   4365.09703163  5394.94930993]
 >>> print dpde
 [ 13441.9  13191.9  13691.9]
"""
import numpy as np
import os
from collections import namedtuple

from xfv.src.equationsofstate.equationofstatebase import EquationOfStateBase

# Deactivate pylint warnings due to NotImplementedError
# pylint: disable=R0921

MieGruneisenParameters = namedtuple('MieGruneisenParameters', ['czero', 'S1', 'S2', 'S3',
                                                               'rhozero', 'grunzero', 'b', 'ezero'])


class MieGruneisen(EquationOfStateBase):
    """
    Mie Gruneisen equation of state
    """

    def __init__(self, czero=3940.0, S1=1.489, S2=0., S3=0., rhozero=8930.0, grunzero=2.02,
                 b=0.47, ezero=0.):
        #
        self.__param = MieGruneisenParameters(czero=czero, S1=S1, S2=S2, S3=S3, rhozero=rhozero,
                                              grunzero=grunzero, b=b, ezero=ezero)
        self.__czero2 = self.__param.czero ** 2
        self.__dgam = self.__param.rhozero * (self.__param.grunzero - self.__param.b)

    def __str__(self):
        message = "EquationOfState : {:s}".format(self.__class__.__name__) + os.linesep
        message += "Parameters : "
        for key, val in zip(self.__param._fields, self.__param):
            message += os.linesep + " -- {:>20s} : {:>9.8g}".format(key, val)
        return message

    def __repr__(self):
        msg = ", ".join(["{:s}={:f}".format(k, v) for k, v in zip(self.__param._fields,
                                                                  self.__param)])
        return "MieGruneisen({:s})".format(msg)

    @property
    def eos_param(self):
        """
        Accesseur sur les param�tres de l'�quation d'�tat
        :return:
        """
        return self.__param

    def solveVolumeEnergy(self, specific_volume, internal_energy, pressure, vson, gampervol):
        """
        Given the specific volume and internal energy computes the pressure, sound speed and
        derivative of the pressure with respect to the internal energy

        :param specific_volume: specific volume (in)
        :type specific_volume: numpy.array
        :param internal_energy: internal energy (in)
        :type internal_energy: numpy.array
        :param pressure: pressure (out)
        :type pressure: numpy.array
        :param vson: sound speed (out)
        :type vson: numpy.array
        :param gampervol: derivative of pressure with respect to the internal energy (out)
        :type gampervol: numpy.array
        """
        epsv = 1 - self.__param.rhozero * specific_volume
        gampervol[:] = 1. / specific_volume
        gampervol *= (self.__param.grunzero * (1 - epsv) + self.__param.b * epsv)
        #
        targets = epsv > 0  # �Cells in compression (~targets are cells in release)
        self.__compression_case(specific_volume, internal_energy, pressure, vson, gampervol, epsv,
                                targets)
        self.__release_case(specific_volume, internal_energy, pressure, vson, gampervol, epsv,
                            ~targets)
        pb = vson < 0
        if pb.any():
            msg = "Sound speed square < 0 in cells {}\n".format(np.where(pb))
            msg += "specific_volume = {}\n".format(specific_volume[pb])
            msg += "pressure = {}\n".format(pressure[pb])
            msg += "dpsurde = {}\n".format(gampervol[pb])
            raise ValueError(msg)
        vson[:] = np.sqrt(vson[:])
        #
        pb = vson >= 10000.
        vson[pb] = 0.

    def __release_case(self, specific_volume, internal_energy, pressure, vson2, gampervol, epsv,
                       targets):
        """
        Compute the equation of state for cells under release conditions.

        :param specific_volume: specific volume (in)
        :type specific_volume: numpy.array
        :param internal_energy: internal energy (in)
        :type internal_energy: numpy.array
        :param pressure: pressure (out)
        :type pressure: numpy.array
        :param vson2: sound speed square (out)
        :type vson2: numpy.array
        :param gampervol: derivative of pressure with respect to the internal energy (out)
        :type gampervol: numpy.array
        :param epsv: compression of the material
        :type epsv: numpy.array
        :param targets: Mask representing the cells in release state
        :type targets: numpy.array (bool)
        """
        loc_epsv = epsv[targets]
        loc_gampervol = gampervol[targets]
        #
        phi = self.__param.rhozero * self.__czero2 * loc_epsv / (1. - loc_epsv)
        # einth ---> e0
        einth = self.__param.ezero
        #
        dphi = -self.__czero2 / specific_volume[targets] ** 2
        #
        dpdv = dphi + (self.__dgam - loc_gampervol) * \
                      (internal_energy[targets] - einth) / specific_volume[targets]
        pressure[targets] = phi + loc_gampervol * (internal_energy[targets] - einth)

        vson2[targets] = specific_volume[targets] ** 2 * (pressure[targets] * loc_gampervol - dpdv)

    def __compression_case(self, specific_volume, internal_energy, pressure, vson2, gampervol,
                           epsv, targets):
        """
        Compute the equation of state for cells under compressive conditions.

        :param specific_volume: specific volume (in)
        :type specific_volume: numpy.array
        :param internal_energy: internal energy (in)
        :type internal_energy: numpy.array
        :param pressure: pressure (out)
        :type pressure: numpy.array
        :param vson2: sound speed square (out)
        :type vson2: numpy.array
        :param gampervol: derivative of pressure with respect to the internal energy (out)
        :type gampervol: numpy.array
        :param epsv: compression of the material
        :type epsv: numpy.array
        :param targets: Mask representing the cells in compressive state
        :type targets: numpy.array (bool)
        """
        loc_epsv = epsv[targets]
        loc_epsv2 = loc_epsv ** 2
        # Coefficient de gruneisen
        loc_gampervol = gampervol[targets]
        redond_a = self.__param.S1 + 2. * self.__param.S2 * loc_epsv \
                   + 3. * self.__param.S3 * loc_epsv2
        denom = (1. - (self.__param.S1 + self.__param.S2 * loc_epsv
                       + self.__param.S3 * loc_epsv2) * loc_epsv)
        phi = self.__param.rhozero * self.__czero2 * loc_epsv / denom ** 2
        einth = self.__param.ezero + phi * loc_epsv / (2. * self.__param.rhozero)
        #
        dphi = phi * self.__param.rhozero * (-1. / loc_epsv - 2. * redond_a / denom)
        #
        deinth = phi * (-1. - loc_epsv * redond_a / denom)
        #
        dpdv = dphi + (self.__dgam - loc_gampervol) * \
                      (internal_energy[targets] - einth) / \
               specific_volume[targets] - loc_gampervol * deinth
        pressure[targets] = phi + loc_gampervol * (internal_energy[targets] - einth)
        vson2[targets] = specific_volume[targets] ** 2 * (pressure[targets] * loc_gampervol - dpdv)

    def solveVolumePressure(self, specific_volume, pressure):
        raise NotImplementedError

    def solveVolumeTemperature(self, specific_volume, temperatue):
        raise NotImplementedError
