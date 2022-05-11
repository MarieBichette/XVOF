# -*- coding: utf-8 -*-
"""
Interface for the porosity calculation
"""
import numpy as np
from xfv.src.porosity_model.porositymodel_base import PorosityModelBase

class JohnsonModel(PorosityModelBase):  # pylint: disable=too-few-public-methods
    """
    Abstract class for the porosity model computation
    """
    def __init__(self,
                 initial_porosity_for_johnson,
                 effective_strength_for_johnson,
                 viscosity_for_johnson, maximal_porosity_for_Johnson):
        self.initial_porosity_for_johnson = initial_porosity_for_johnson
        self.effective_strength_for_johnson = effective_strength_for_johnson
        self.viscosity_for_johnson = viscosity_for_johnson
        self.maximal_porosity_for_Johnson = maximal_porosity_for_Johnson

    def compute_porosity(self, delta_t: float,
                         porosity: np.array,
                         pressure: np.array)->np.array:
        """
        Check that the tension is greater than the equilibium tension and
        compute the new value of porosity
        """
        initial_porosity_for_johnson = self.initial_porosity_for_johnson
        effective_strength_for_johnson = self.effective_strength_for_johnson

        peq0 = np.ones_like(pressure)*(effective_strength_for_johnson /
                                       initial_porosity_for_johnson *
                                       np.log(initial_porosity_for_johnson /
                                              (initial_porosity_for_johnson-1.0)))

        is_pressure_above_peq = -pressure >= peq0
        is_porosity_one = porosity == 1.0
        is_poro_evolution_ok = np.logical_and(is_pressure_above_peq, is_porosity_one)
        porosity[is_poro_evolution_ok] = initial_porosity_for_johnson

        tmp_delta = (abs(pressure)-effective_strength_for_johnson /
                     porosity * np.log(porosity/(porosity-1.0)))
        delta_p = np.where(porosity > 1.0, tmp_delta, 0.)
        delta_p = np.maximum(delta_p, 0.)

        return self._compute_johnson_porosity(delta_t, pressure, porosity, delta_p)

    def _compute_johnson_porosity(self, delta_t: float,
                                  pressure: np.array,
                                  porosity: np.array,
                                  delta_p: np.array)->np.array:
        """
        Compute the new value of porosity
        """
        initial_porosity_for_johnson = self.initial_porosity_for_johnson
        maximal_porosity_for_Johnson = self.maximal_porosity_for_Johnson
        viscosity_for_johnson = self.viscosity_for_johnson

        power = 1.0/3.0
        pexp = np.power(np.array(porosity-1.0), power)

        dalphadt = -((initial_porosity_for_johnson - 1.0)**(2.0/3.0)*
                     porosity * pexp*delta_p) * np.sign(pressure)
        dalphadt = dalphadt/viscosity_for_johnson
        porosity_new = porosity + dalphadt*delta_t
        porosity_new = np.maximum(porosity_new, initial_porosity_for_johnson)
        porosity_new = np.minimum(porosity_new, maximal_porosity_for_Johnson)
        porosity_new = np.where(porosity_new > 1., porosity_new, 1.)
        return porosity_new
