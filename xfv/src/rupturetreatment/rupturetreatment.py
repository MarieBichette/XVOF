# -*- coding: iso-8859-1 -*-
"""
Defining an interface for rupture treatments
"""
from abc import ABCMeta, abstractmethod


class RuptureTreatment(object, metaclass=ABCMeta):  # pylint: disable=too-few-public-methods
    """
    An interface for rupture treatments
    """

    @abstractmethod
    def apply_treatment(self, cells, ruptured_cells, *args, **kwargs):
        """
        Application of the treatment on the ruptured cells

        :param cells: array of all cells
        :param ruptured_cells: boolean array marking the ruptured cells
        """
        pass
