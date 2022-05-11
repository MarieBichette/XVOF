#!/usr/bin/env python3.7
# -*- coding: iso-8859-1 -*-
"""
Implementing a non local criterion for failure + the cell is cracked if it verifies the failure criterion itself
"""
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion
from xfv.src.rupturecriterion.nonlocalstress import NonLocalStressCriterion


class NonLocalStressCriterionWithMax(RuptureCriterion):  # pylint: disable=too-few-public-methods
    """
    An class for non local failure criterion
    """
    def __init__(self, value, average_strategy):
        super().__init__()
        self.critical_value = value
        self.non_local_criterion = NonLocalStressCriterion(value, average_strategy)

    def check_criterion(self, cells, *args, **kwargs):
        """
        Check of the rupture criterion on the cells in arguments
        :param cells: cells on which to check the criterion
        """
        return self.non_local_criterion.check_criterion(cells) * (cells.stress[:, 0] >= self.critical_value)
