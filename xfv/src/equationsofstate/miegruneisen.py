#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
 >>> my_eos.solve_volume_energy(specific_volume, internal_energy, pressure, dpde, sound_speed)
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

    def solve_volume_energy(self, specific_volume, internal_energy, pressure,  # pylint: disable=too-many-arguments
                            derivative, vson=None):
        """
        Given the specific volume and internal energy computes the pressure, sound speed and
        derivative of the pressure with respect to the internal energy

        :param specific_volume: specific volume (in)
        :type specific_volume: numpy.array
        :param internal_energy: internal energy (in)
        :type internal_energy: numpy.array
        :param pressure: pressure (out)
        :type pressure: numpy.array
        :param derivative: derivative of pressure with respect to the internal energy (out)
        :type derivative: numpy.array
        :param vson: sound speed (out)
        :type vson: numpy.array
        """
        epsv = 1 - self.__param.rhozero * specific_volume
        derivative[:] = 1. / specific_volume
        derivative *= (self.__param.grunzero * (1 - epsv) + self.__param.b * epsv)
        einth = np.ndarray(specific_volume.shape)
        phi = np.ndarray(specific_volume.shape)
        #
        comp_cells = epsv > 0  # Cells in compression
        rel_cells = ~comp_cells  # Cells in release
        epsv_comp = epsv[comp_cells]
        epsv_rel = epsv[rel_cells]
        denom = np.ndarray(epsv_comp.shape)
        einth[comp_cells], phi[comp_cells] = self.__compute_eint_phi_compression(epsv_comp, denom)
        einth[rel_cells], phi[rel_cells] = self.__compute_eint_phi_release(epsv_rel)
        pressure[:] = phi + derivative * (internal_energy - einth)

        if vson is not None:
            self.__compute_vson(specific_volume, internal_energy, pressure, derivative,
                                epsv, einth, phi, denom, comp_cells, rel_cells, vson)

    def __compute_denom_1(self, epsv):
        return 1. - self.__param.S1 * epsv

    def __compute_denom_2(self, epsv):
        return self.__compute_denom_1(epsv) - self.__param.S2 * epsv * epsv

    def __compute_denom_3(self, epsv):
        return self.__compute_denom_2(epsv) - self.__param.S3 * epsv * epsv * epsv

    def __compute_redonda_1(self):
        return self.__param.S1

    def __compute_redonda_2(self, epsv):
        return self.__compute_redonda_1() + 2. * self.__param.S2 * epsv

    def __compute_redonda_3(self, epsv):
        return self.__compute_redonda_2(epsv) + 3. * self.__param.S3 * epsv ** 2

    def __compute_eint_phi_compression(self, epsv, denom):
        if self.__param.S2 and self.__param.S3:
            denom[:] = self.__compute_denom_3(epsv)
        elif self.__param.S2:
            denom[:] = self.__compute_denom_2(epsv)
        else:
            denom[:] = self.__compute_denom_1(epsv)
        phi = self.__param.rhozero * self.__czero2 * epsv / denom ** 2
        einth = self.__param.ezero + phi * epsv / (2. * self.__param.rhozero)
        return einth, phi

    def __compute_eint_phi_release(self, epsv):
        phi = self.__param.rhozero * self.__czero2 * epsv / (1. - epsv)
        einth = self.__param.ezero
        return einth, phi

    def __compute_dpdv_compression(self, specific_volume, internal_energy,  # pylint: disable=too-many-arguments, too-many-locals
                                   gampervol, epsv, einth, phi, denom):
        if self.__param.S2 and self.__param.S3:
            redond_a = self.__compute_redonda_3(epsv)
        elif self.__param.S2:
            redond_a = self.__compute_redonda_2(epsv)
        else:
            redond_a = self.__compute_redonda_1()
        dphi = phi * self.__param.rhozero * (-1. / epsv - 2. * redond_a / denom)
        deinth = phi * (-1. - epsv * redond_a / denom)
        dpdv = (dphi + (self.__dgam - gampervol) *
                (internal_energy - einth) / specific_volume - gampervol * deinth)
        return dpdv

    def __compute_dpdv_release(self, specific_volume, internal_energy,  # pylint: disable=too-many-arguments
                               gampervol, einth):
        dphi = - self.__czero2 / specific_volume ** 2
        dpdv = (dphi + (self.__dgam - gampervol) *
                (internal_energy - einth) / specific_volume)
        return dpdv

    def __compute_vson(self, specific_volume, internal_energy, pressure, derivative,  # pylint: disable=too-many-arguments, too-many-locals
                       epsv, einth, phi, denom, comp_cells, rel_cells, vson):
        dpdv = np.ndarray(specific_volume.shape)
        specific_volume_comp = specific_volume[comp_cells]
        specific_volume_rel = specific_volume[rel_cells]
        internal_energy_comp = internal_energy[comp_cells]
        internal_energy_rel = internal_energy[rel_cells]
        derivative_comp = derivative[comp_cells]
        derivative_rel = derivative[rel_cells]
        epsv_comp = epsv[comp_cells]
        einth_comp = einth[comp_cells]
        einth_rel = einth[rel_cells]
        phi_comp = phi[comp_cells]
        dpdv[comp_cells] = self.__compute_dpdv_compression(
            specific_volume_comp, internal_energy_comp, derivative_comp,
            epsv_comp, einth_comp, phi_comp, denom)
        dpdv[rel_cells] = self.__compute_dpdv_release(
            specific_volume_rel, internal_energy_rel, derivative_rel, einth_rel)
        vson[:] = specific_volume ** 2 * (pressure * derivative - dpdv)
        vson[:] = np.sqrt(vson)
        vson[vson >= 10000.] = 0.

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
