#!/usr/bin/env python3.7
# -*- coding: iso-8859-1 -*-
"""
Implementing a computation of the weight associated to each neighbour cell
"""
import numpy as np
from abc import ABCMeta, abstractmethod


class IWeight(metaclass=ABCMeta):
    def __init__(self, radius):
        self.radius = radius
        pass

    @abstractmethod
    def compute_weight(self, length):
        return 0


class ArithmeticWeight(IWeight):
    def ___init__(self, radius):
        """
        A class for no weight (weight = 1 for each neighbour)
        """
        super(ArithmeticWeight, self).__init__(radius)

    def compute_weight(self, length):
        """
        Weight = 1 for cells in the neighbourhood
        """
        weight = np.zeros(length.shape)
        mask_in = np.abs(length) <= self.radius
        weight[mask_in] = len(np.where(mask_in)[0]) * [1.]
        return weight


class LinearWeight(IWeight):
    def __init__(self, radius):
        """
        A class for linear weight
        """
        super(LinearWeight, self).__init__(radius)

    def compute_weight(self, length):
        """
        The weight is linearly decreasing with the distance to cell
        """
        weight = np.zeros(length.shape)
        mask_in = np.abs(length) <= self.radius
        weight[mask_in] = 1. - np.abs(length[mask_in] / self.radius)
        return weight


class GaussianWeight(IWeight):
    def __init__(self, radius):
        """
        A class for gaussian weight
        """
        super(GaussianWeight, self).__init__(radius)

    def compute_weight(self, length):
        """
        The weight is exponentially decreasing with the distance to cell
        """
        raise NotImplementedError("""Il manque des param�tres dans lejeu de donn�es pour calculer la gaussienne""")

