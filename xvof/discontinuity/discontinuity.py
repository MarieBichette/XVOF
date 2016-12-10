# -*- coding: iso-8859-1 -*-
"""
A module implementing the Discontinuity class
"""


class Discontinuity(object):
    """
    A class describing a discontinuity
    """

    discontinuity_number = 0

    def __init__(self, mask_in_nodes, mask_out_nodes):
        Discontinuity.discontinuity_number += 1
        self.__label = Discontinuity.discontinuity_number
        self.__mask_in_nodes = mask_in_nodes
        self.__mask_out_nodes = mask_out_nodes
        print "Building discontinuity number {:d}".format(self.__label)

    @property
    def mask_in_nodes(self):
        """
        Accessor on the mask on the nodes "in" the discontinuity

        :return: the mask on the nodes "in" the discontinuity
        """
        return self.__mask_in_nodes

    @property
    def mask_out_nodes(self):
        """
        Accessor on the mask on the nodes "out" the discontinuity

        :return: the mask on the nodes "out" the discontinuity
        """
        return self.__mask_out_nodes
