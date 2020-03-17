#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module pour calculer les invariants du tenseur des contraintes
"""
import numpy as np


def compute_trace(stress):
    """
    Compute trace(sigma)
    :return : trace(sigma)
    """

    trace = stress[:, 0] + stress[:, 1] + stress[:, 2]
    return trace


def compute_J2(stress):
    """
    Compute the second invariant of stress tensor
    J2 = sqrt(3/2 S:S)
    where S is the deviatoric part of stress tensor
    """
    trace = compute_trace(stress)
    # on calcule de déviateur de stress (pas d'effet si stress est deviatoric_stress_...)
    S = np.copy(stress)
    S[:, 0] -= 1. / 3. * trace
    S[:, 1] -= 1. / 3. * trace
    S[:, 2] -= 1. / 3. * trace
    # en 1D : pour l'instant le tenseur S est diagonal
    # Il est enregistré comme un vecteur de taille 3 pour chaque élément = array de taille nbr of element * 3
    # le calcul de J2 prend en compte cette caractéristique (pas de terme de cisaillement)
    J2 = S[:, 0]**2. + S[:, 1]**2 + S[:, 2]**2
    J2 *= 3./2.
    J2 = np.sqrt(J2)
    return J2