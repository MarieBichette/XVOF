#!/usr/bin/env python3.7
# -*- coding: iso-8859-1 -*-
"""
Implementing a computation of the weight associated to each neighbour cell
"""
import numpy as np
from abc import ABCMeta, abstractmethod


class IWeight(metaclass=ABCMeta):
    def __init__(self):
        pass

    @abstractmethod
    def compute_weight(self, radius, length):
        return 0


class NoWeight(IWeight):
    def ___init__(self):
        """
        A class for no weight (weight = 1 for each neighbour)
        """
        super(NoWeight, self).__init__()

    def compute_weight(self, radius, length):
        """
        Weight = 1
        """
        weight = np.zeros_like(length)
        mask_in = np.abs(length) > radius
        weight[mask_in] = np.ones_like(length[mask_in])
        return weight


class LinearWeight(IWeight):
    def __init__(self):
        """
        A class for linear weight
        """
        super(LinearWeight, self).__init__()

    def compute_weight(self, radius, length):
        """
        The weight is linearly decreasing with the distance to cell
        """
        weight = np.zeros_like(length)
        mask_in = np.abs(length) > radius
        weight[mask_in] = 1. - np.abs(length[mask_in] / radius)
        return weight


class GaussianWeight(IWeight):
    def __init__(self):
        """
        A class for gaussian weight
        """
        super(GaussianWeight, self).__init__()

    def compute_weight(self, radius, length):
        """
        The weight is exponentially decreasing with the distance to cell
        """
        raise NotImplementedError("""Il manque des paramètres dans lejeu de données pour calculer la gaussienne""")

