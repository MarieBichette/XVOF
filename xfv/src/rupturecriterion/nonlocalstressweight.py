#!/usr/bin/env python3.7
# -*- coding: iso-8859-1 -*-
"""
Implementing a computation of the weight associated to each neighbour cell
"""
import numpy as np
from abc import abstractmethod


class IWeight:
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
        return 1


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
        if np.abs(length) > radius:
            return 0
        return 1. - np.abs(length) / radius


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
        if np.abs(length) > radius:
            return 0
        return 1. - np.exp(np.abs(length) / radius)
