#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module OneDimensionEnrichedNode
"""
import unittest
import mock
import os
import numpy as np

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
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_hydro.xml")
        self.test_datacontainer = DataContainer(data_file_path)

        self.vit_init = np.zeros([4, 1], dtype='float')
        self.vit_init[:, 0] = [-1.5e+03, 1.2e+03, 0.0, 0.3e+03]
        self.poz_init = np.zeros([4, 1], dtype='float')
        self.poz_init[:, 0] = [0., 1., 2., 3.]
        self.my_nodes = OneDimensionEnrichedNode(4, self.poz_init, self.vit_init, section=1.0e-06)
        self.my_nodes._classical = np.array([True, False, False, True])

        # configuration propre d'un mock 'discontinuity'
        config = {'mask_in_nodes': np.array([False, True, False, False]),
                  'mask_out_nodes': np.array([False, False, True, False]),
                  'additional_dof_velocity_current': np.array([[1., ], [1., ]]),
                  'additional_dof_force': np.array([[1., ], [2., ]]),
                  'label': 1
                  }
        patcher = mock.patch('xfv.src.discontinuity.discontinuity.Discontinuity',
                             spec=Discontinuity, **config)
        self.mock_discontinuity = patcher.start()

    def tearDown(self):
        """
        Operations to be done after completing all the tests in the class
        """
        DataContainer.clear()

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_additional_dof_new_velocity(self, mock_disc_list):
        Discontinuity.discontinuity_list.return_value = [self.mock_discontinuity]
        inv_mass_additional = np.array([[2., 0.], [0., 2.]])
        self.my_nodes.compute_additional_dof_new_velocity(1., inv_mass_additional)
        np.testing.assert_array_almost_equal(self.mock_discontinuity._additional_dof_velocity_new,
                                             np.array([[3., ], [5., ]]))


if __name__ == '__main__':
    unittest.main()
