""" Test of wc_test.core

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-05-10
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_test
import unittest

class TestCore(unittest.TestCase):

    def test_1(self):
        model_path = '/home/balazs/Desktop/wc_test/tests/fixtures/min_model.xlsx'
        self.assertFalse(wc_test.SubmodelDynamicsTestCase().is_species_constant(model_path, 600, 10, 'ATP', 'c', 'RNA_1', 'c'))
        self.assertFalse(wc_test.SubmodelDynamicsTestCase().is_species_constant(model_path, 600, 10, 'ATP', 'c', 'RNA_2', 'c'))
        self.assertFalse(wc_test.SubmodelDynamicsTestCase().is_species_constant(model_path, 600, 10, 'ATP', 'c', 'RNA_3', 'c'))
        self.assertTrue(wc_test.SubmodelDynamicsTestCase().is_species_constant(model_path, 600, 10, 'ATP', 'c', 'RNA_4', 'c'))
        self.assertTrue(wc_test.SubmodelDynamicsTestCase().is_species_constant(model_path, 600, 10, 'ATP', 'c', 'RNA_5', 'c'))
