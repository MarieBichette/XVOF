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
 MieGruneisen(czero=3980.000000, S1=1.580000, S2=0.000000, S3=0.000000,
 rhozero=8129.000000, grunzero=1.600000, b=0.500000, ezero=0.000000)
 >>> density = np.array([9000., 8500., 9500.], dtype=np.float64, order='C')
 >>> specific_volume = 1. / density
 >>> internal_energy = np.array([1.0e+04, 1.0e+03, 1.0e+05], dtype=np.float64, order='C')
 >>> the_shape = internal_energy.shape
 >>> pressure = np.zeros(the_shape, dtype=np.float64, order='C')
 >>> sound_speed = np.zeros(the_shape, dtype=np.float64, order='C')
 >>> dpde = np.zeros(the_shape, dtype=np.float64, order='C')
 >>> my_eos.solve_volume_energy(specific_volume, internal_energy, pressure, sound_speed, dpde)
 >>> print pressure
 [  1.61115797e+10   6.26727977e+09   2.87613980e+10]
 >>> print sound_speed
 [ 4871.9323597   4365.09703163  5394.94930993]
 >>> print dpde
 [ 13441.9  13191.9  13691.9]
"""
import os
from collections import namedtuple

import numpy as np

from xfv.src.equationsofstate.equationofstatebase import EquationOfStateBase

# Deactivate pylint warnings due to NotImplementedError

MieGruneisenParameters = namedtuple('MieGruneisenParameters', ['czero', 'S1', 'S2', 'S3',
                                                               'rhozero', 'grunzero', 'b', 'ezero'])


class MieGruneisen(EquationOfStateBase):
    """
    Mie Gruneisen equation of state
    """

    def __init__(self, czero=3940.0, S1=1.489, S2=0., S3=0., rhozero=8930.0,  # pylint: disable=too-many-arguments
                 grunzero=2.02, b=0.47, ezero=0.):
        #
        self.__param = MieGruneisenParameters(czero=czero, S1=S1, S2=S2, S3=S3,
                                              rhozero=rhozero, grunzero=grunzero,
                                              b=b, ezero=ezero)
        self.__czero2 = self.__param.czero ** 2
        self.__dgam = self.__param.rhozero * (self.__param.grunzero - self.__param.b)

    def __str__(self):
        message = "EquationOfState : {:s}".format(self.__class__.__name__) + os.linesep
        message += "Parameters : "
        for key, val in zip(self.__param._fields, self.__param):
            message += os.linesep + " -- {:>20s} : {:>9.8g}".format(key, val)
        return message

    def __repr__(self):
        msg = ", ".join(["{:s}={:f}".format(k, v)
                         for k, v in zip(self.__param._fields, self.__param)])
        return "MieGruneisen({:s})".format(msg)

    @property
    def eos_param(self):
        """
        Accesseur sur les param�tres de l'�quation d'�tat
        :return:
        """
        return self.__param

    def solve_volume_energy(self, specific_volume, internal_energy, pressure, vson, derivative):  # pylint: disable=too-many-arguments
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
        derivative[:] = 1. / specific_volume
        derivative *= (self.__param.grunzero * (1 - epsv) + self.__param.b * epsv)
        dpdv = np.ndarray(specific_volume.shape)
        #
        targets = epsv > 0  # �Cells in compression (~targets are cells in release)
        self.__compression_case(specific_volume, internal_energy, pressure, dpdv,
                                derivative, epsv, targets)
        self.__release_case(specific_volume, internal_energy, pressure, dpdv,
                            derivative, epsv, ~targets)
        self.__compute_vson2(specific_volume, pressure, derivative, dpdv, vson)

    def __release_case(self, specific_volume, internal_energy, pressure, dpdv,  # pylint: disable=too-many-arguments
                       gampervol, epsv, targets):
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
        dpdv[targets] = (dphi + (self.__dgam - loc_gampervol) *
                         (internal_energy[targets] - einth) / specific_volume[targets])
        pressure[targets] = phi + loc_gampervol * (internal_energy[targets] - einth)

    @staticmethod
    def __compute_vson2(specific_volume, pressure, gampervol, dpdv, vson):
        """
        Compute the sound speed velocity
        """
        vson[:] = specific_volume[:] ** 2 * (pressure[:] * gampervol[:] - dpdv[:])
        vson[:] = np.sqrt(vson[:])
        vson[vson >= 10000.] = 0.

    def __compression_case(self, specific_volume, internal_energy, pressure, dpdv,  # pylint: disable=too-many-arguments, too-many-locals
                           gampervol, epsv, targets):
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
        redond_a = (self.__param.S1 + 2. * self.__param.S2 * loc_epsv +
                    3. * self.__param.S3 * loc_epsv2)
        denom = (1. -
                 (self.__param.S1 + self.__param.S2 * loc_epsv + self.__param.S3 * loc_epsv2)
                 * loc_epsv)
        phi = self.__param.rhozero * self.__czero2 * loc_epsv / denom ** 2
        einth = self.__param.ezero + phi * loc_epsv / (2. * self.__param.rhozero)
        #
        dphi = phi * self.__param.rhozero * (-1. / loc_epsv - 2. * redond_a / denom)
        #
        deinth = phi * (-1. - loc_epsv * redond_a / denom)
        #
        dpdv[targets] = (dphi +
                         (self.__dgam - loc_gampervol) *
                         (internal_energy[targets] - einth) / specific_volume[targets]
                         - loc_gampervol * deinth)


if __name__ == "__main__":
    import time
    # pylint: disable=invalid-name
    eos = MieGruneisen()
    density = np.ndarray((1000,))
    int_nrj = np.ndarray((1000,))
    press = np.ndarray((1000,))
    sound_speed = np.ndarray((1000,))
    deriv = np.ndarray((1000,))

    density[:] = 9000.
    int_nrj[:] = 1.e+05
    press[:] = 0.
    sound_speed[:] = 0.
    deriv[:] = 0.

    spec_vol = 1. / density
    start_time = time.perf_counter()
    iter_number = 15000
    for _ in range(iter_number):
        eos.solve_volume_energy(spec_vol, int_nrj, press, deriv, sound_speed)
    end_time = time.perf_counter()

    print(f"{iter_number} iterations : {end_time - start_time} seconds")
