""" Test of wc_test.core

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-05-10
:Copyright: 2018, Karr Lab
:License: MIT

Todo: add more cases to make sure every  edge case is covered
"""

import wc_lang
import wc_test
import unittest

class TestCore(unittest.TestCase):
    model_path = 'tests/fixtures/min_model.xlsx'
    target_specie_ids = ['RNA_1[c]','RNA_2[c]','RNA_3[c]','RNA_4[c]','RNA_5[c]']

    def test_is_constant_tweak_species(self):
        tweak_specie_ids = ['ATP[c]']
        is_constant, run_results = wc_test.SubmodelDynamicsTestCase().is_constant_tweak_species(model_path=self.model_path,
                  end_time=600, checkpoint_period=10, tweak_specie_ids=tweak_specie_ids, target_specie_ids=self.target_specie_ids)
        self.assertEqual(is_constant, ['False', 'False', 'False', 'True', 'True'])
        wc_lang.SpeciesType.objects.reset() # reset indexer

        tweak_specie_ids = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']
        is_constant, run_results = wc_test.SubmodelDynamicsTestCase().is_constant_tweak_species(model_path=self.model_path,
                  end_time=600, checkpoint_period=10, tweak_specie_ids=tweak_specie_ids, target_specie_ids=self.target_specie_ids)
        self.assertEqual(is_constant, ['False', 'False', 'False', 'True', 'True']) # Degradation active for first 3 species
        wc_lang.SpeciesType.objects.reset() # reset indexer

    def test_is_constant_tweak_reaction(self):
        tweak_reaction_ids = ['transcription_RNA_4']
        target_specie_ids = ['RNA_1[c]','RNA_2[c]','RNA_3[c]','RNA_4[c]','RNA_5[c]']
        is_constant, run_results = wc_test.SubmodelDynamicsTestCase().is_constant_tweak_reactions(model_path=self.model_path, end_time=600, checkpoint_period=10,
                  tweak_reaction_ids=tweak_reaction_ids, target_specie_ids=self.target_specie_ids)
        self.assertEqual(is_constant, ['False', 'False', 'False', 'True', 'True'])
        wc_lang.SpeciesType.objects.reset() # reset indexer

        tweak_reaction_ids = ['transcription_RNA_1', 'transcription_RNA_2', 'transcription_RNA_3', 'transcription_RNA_4', 'transcription_RNA_5']
        is_constant, run_results = wc_test.SubmodelDynamicsTestCase().is_constant_tweak_reactions(model_path=self.model_path, end_time=600, checkpoint_period=10,
                  tweak_reaction_ids=tweak_reaction_ids, target_specie_ids=self.target_specie_ids)
        self.assertEqual(is_constant, ['False', 'False', 'False', 'True', 'True']) # Degradation active for first 3 species
        wc_lang.SpeciesType.objects.reset() # reset indexer
