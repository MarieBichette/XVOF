#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Module to compute the second invariant of the stress tensor
"""
import numpy as np


def compute_trace(stress: np.array):
    """
    Compute trace(sigma)
    :param stress : stress tensor
    :return : trace(sigma)
    """
    trace = stress[:, 0] + stress[:, 1] + stress[:, 2]
    return trace


def compute_second_invariant(stress: np.array):
    """
    Compute the second invariant of stress tensor
    J2 = sqrt(3/2 S:S) where S is the deviatoric part of stress tensor
    :param stress : stress tensor
    """
    trace = compute_trace(stress)
    # Compute the deviatoric stress
    dev_stress = np.copy(stress)
    dev_stress[:, 0] -= 1. / 3. * trace
    dev_stress[:, 1] -= 1. / 3. * trace
    dev_stress[:, 2] -= 1. / 3. * trace
    second_invariant = dev_stress[:, 0]**2. + dev_stress[:, 1]**2 + dev_stress[:, 2]**2
    second_invariant *= 3./2.
    second_invariant = np.sqrt(second_invariant)
    return second_invariant
