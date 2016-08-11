# -*- coding: iso-8859-1 -*-
"""
Implementing MieGruneisen class

>>> import numpy as np
>>> from xvof.equationsofstate.miegruneisen import MieGruneisen
>>> my_eos = MieGruneisen()
>>> print my_eos # doctest:+NORMALIZE_WHITESPACE
EquationOfState : MieGruneisen
Parameters :
 --                   S1 :      1.58
 --                   S2 :         0
 --                   S3 :         0
 --                    b :       0.5
 --                czero :      3980
 --               czero2 :  15840400
 --                 dgam :    8941.9
 --                ezero :         0
 --             grunzero :       1.6
 --              rhozero :      8129
 >>> my_eos # doctest:+NORMALIZE_WHITESPACE
 MieGruneisen(S1=1.580000, S2=0.000000, S3=0.000000, b=0.500000, czero=3980.000000, czero2=15840400.000000,
 dgam=8941.900000, ezero=0.000000, grunzero=1.600000, rhozero=8129.000000)
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
from UserDict import UserDict
from operator import getitem

from xvof.equationsofstate.equationofstatebase import EquationOfStateBase


# Deactivate pylint warnings due to NotImplementedError
# pylint: disable=R0921


class MieGruneisen(UserDict, EquationOfStateBase):
    """
    Mie Gruneisen equation of state
    """

    def __init__(self, czero=3980.0, S1=1.58, S2=0., S3=0., rhozero=8129.0, grunzero=1.6, b=0.5, ezero=0.):
        #
        super(MieGruneisen, self).__init__(czero=czero, S1=S1, S2=S2, S3=S3, rhozero=rhozero, grunzero=grunzero, b=b,
                                           ezero=ezero, czero2=czero ** 2, dgam=rhozero * (grunzero - b))

    def __str__(self):
        message = "EquationOfState : {:s}".format(self.__class__.__name__)
        message += "\nParameters : "
        for key, val in sorted(self.iteritems(), key=lambda m: getitem(m, 0)):
            message += "\n -- {:>20s} : {:>9.8g}".format(key, val)
        return message

    def __repr__(self):
        msg = ", ".join(["{:s}={:f}".format(k, v) for k, v in sorted(self.iteritems(), key=lambda m: getitem(m, 0))])
        return "MieGruneisen({:s})".format(msg)

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
        epsv = 1 - self["rhozero"] * specific_volume
        gampervol[:] = (self["grunzero"] * (1 - epsv) + self["b"] * epsv) / specific_volume
        #
        targets = epsv > 0  #  Cells in compression (~targets are cells in release)
        self.__compression_case(specific_volume, internal_energy, pressure, vson, gampervol, epsv, targets)
        self.__release_case(specific_volume, internal_energy, pressure, vson, gampervol, epsv, ~targets)
        vson[:] = np.sqrt(vson)
        #
        pb = vson >= 10000.
        vson[pb] = 0.

    def __release_case(self, specific_volume, internal_energy, pressure, vson2, gampervol, epsv, targets):
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
        phi = self["rhozero"] * self["czero2"] * loc_epsv / (1. - loc_epsv)
        # einth ---> e0
        einth = self["ezero"]
        #
        dphi = -self["czero2"] / specific_volume[targets] ** 2
        #
        dpdv = dphi + (self["dgam"] - loc_gampervol) * \
                      (internal_energy[targets] - einth) / specific_volume[targets]
        pressure[targets] = phi + loc_gampervol * (internal_energy[targets] - einth)
        vson2[targets] = specific_volume[targets] ** 2 * (pressure[targets] * loc_gampervol - dpdv)
        pb = vson2 < 0
        if pb.any():
            msg = "Sound speed square < 0\n"
            msg += "specific_volume = {}\n".format(specific_volume[pb])
            msg += "pressure = {}\n".format(pressure[pb])
            msg += "dpsurde = {}\n".format(gampervol[pb])
            msg += "dpdv = {}\n".format(dpdv)
            raise ValueError(msg)

    def __compression_case(self, specific_volume, internal_energy, pressure, vson2, gampervol, epsv, targets):
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
        redond_a = self["S1"] + 2. * self["S2"] * loc_epsv + 3. * self["S3"] * loc_epsv2
        denom = (1. - (self["S1"] + self["S2"] * loc_epsv + self["S3"] * loc_epsv2) * loc_epsv)
        phi = self["rhozero"] * self["czero2"] * loc_epsv / denom ** 2
        einth = self["ezero"] + phi * loc_epsv / (2. * self["rhozero"])
        #
        dphi = phi * self["rhozero"] * (-1. / loc_epsv - 2. * redond_a / denom)
        #
        deinth = phi * (-1. - loc_epsv * redond_a / denom)
        #
        dpdv = dphi + (self["dgam"] - loc_gampervol) * \
                      (internal_energy[targets] - einth) / specific_volume[targets] - loc_gampervol * deinth
        pressure[targets] = phi + loc_gampervol * (internal_energy[targets] - einth)
        vson2[targets] = specific_volume[targets] ** 2 * (pressure[targets] * loc_gampervol - dpdv)
        pb = vson2 < 0
        if pb.any():
            msg = "Sound speed square < 0\n"
            msg += "specific_volume = {}\n".format(specific_volume[pb])
            msg += "pressure = {}\n".format(pressure[pb])
            msg += "dpsurde = {}\n".format(gampervol[pb])
            msg += "dpdv = {}\n".format(dpdv)
            raise ValueError(msg)

    def solveVolumePressure(self, specific_volume, pressure):
        raise NotImplementedError

    def solveVolumeTemperature(self, specific_volume, temperatue):
        raise NotImplementedError
