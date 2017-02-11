# -*- coding: iso-8859-1 -*-
"""
A module implementing the Discontinuity class
"""


class Discontinuity(object):
    """
    A class describing a discontinuity
    """

    # A list of discontinuities
    # modified in rupturetreatment.enrichelement(...append)
    __discontinuity_list = []

    def __init__(self, mask_in_nodes, mask_out_nodes):
        Discontinuity.__discontinuity_list.append(self)
        self.__label = len(Discontinuity.__discontinuity_list)
        self.__mask_in_nodes = mask_in_nodes
        self.__mask_out_nodes = mask_out_nodes
        self.__mass_matrix_updated = False
        print "Building discontinuity number {:d}".format(self.__label)

    @classmethod
    def discontinuity_number(cls):
        return len(Discontinuity.__discontinuity_list)

    @classmethod
    def discontinuity_list(cls):
        return Discontinuity.__discontinuity_list

    @property
    def position_in_ruptured_element(self):
        """
        Accessor on the "alpha" position in ruptured element
        """
        alpha = 0.5
        return alpha


    @property
    def label(self):
        """
        Accessor on label variable
        :return: label
        """
        return self.__label

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

    @property
    def mass_matrix_updated(self):
        """
        Accessor on the boolean that indicates if the mass matrix has been computed for this discontinuity

        :return: he boolean that indicates if the mass matrix has been computed for this discontinuity
        """
        return self.__mass_matrix_updated

    def hasMassMatrixBeenComputed(self):
        """
        Set the __mass_matrix_updated boolean to True
        """
        self.__mass_matrix_updated = True