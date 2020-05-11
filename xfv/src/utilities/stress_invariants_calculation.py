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


def compute_second_invariant(dev_stress: np.array):
    """
    Compute the square of the second invariant of stress tensor
    J2 = sqrt(3/2 S:S) where S is the deviatoric part of stress tensor
    :param stress : deviatoric stress tensor
    """
    second_invariant = np.multiply(dev_stress, dev_stress)
    second_invariant = second_invariant[:, 0] + second_invariant[:, 1] + second_invariant[:, 2]
    second_invariant *= 3./2.
    return second_invariant
