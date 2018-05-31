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
        tweak_specie_ids = ['ATP[c]']
        target_specie_ids = ['RNA_1[c]','RNA_2[c]','RNA_3[c]','RNA_4[c]','RNA_5[c]']
        results = wc_test.SubmodelDynamicsTestCase().is_species_constant(model_path=model_path, end_time=600, checkpoint_period=10,
                                                               tweak_specie_ids=tweak_specie_ids, target_specie_ids=target_specie_ids)

        self.assertEqual(results, ['False', 'False', 'False', 'True', 'True'])
