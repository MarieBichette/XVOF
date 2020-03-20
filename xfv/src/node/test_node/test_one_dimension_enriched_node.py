#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module OneDimensionEnrichedNode
"""
import numpy as np
import unittest
import mock
import os

from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.node.one_dimension_enriched_node import OneDimensionEnrichedNode
from xfv.src.data.data_container import DataContainer


class OneDimensionEnrichedNodeTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'OneDimensionEnrichedNode'
    """
    def setUp(self):
        """
        Préparation des tests unitaires
        """
        data_file_path = os.path.join(os.path.dirname(__file__), "../../../tests/0_UNITTEST/XDATA_hydro.xml")
        self.test_datacontainer = DataContainer(data_file_path)

        self.vit_init = np.zeros([4, 1], dtype='float')
        self.vit_init[:, 0] = [-1.5e+03, 1.2e+03, 0.0, 0.3e+03]
        self.poz_init = np.zeros([4, 1], dtype='float')
        self.poz_init[:, 0] = [0., 1., 2., 3.]
        self.my_nodes = OneDimensionEnrichedNode(4, self.poz_init, self.vit_init, section=1.0e-06)
        self.my_nodes._classical = np.array([True, False, False, True])

    def tearDown(self):
        DataContainer.clear()
        pass

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_additional_dof_new_velocity(self, mock_disc_list):
        disc = mock.MagicMock(Discontinuity)
        type(disc).mask_in_nodes= mock.PropertyMock(return_value=np.array([False, True, False, False]))
        type(disc).mask_out_nodes = mock.PropertyMock(return_value=np.array([False, False, True, False]))
        type(disc).additional_dof_velocity_current = mock.PropertyMock(return_value=np.array([[1., ], [1., ]]))
        type(disc).additional_dof_force = mock.PropertyMock(return_value=np.array([[1., ], [2., ]]))
        Discontinuity.discontinuity_list.return_value = [disc]

        inv_mass_additional = np.array([[2., 0.],[0., 2.]])
        self.my_nodes.compute_additional_dof_new_velocity(1., inv_mass_additional)
        np.testing.assert_array_almost_equal(disc._additional_dof_velocity_new, np.array([[3., ], [5., ]]))


if __name__ == '__main__':
    unittest.main()
