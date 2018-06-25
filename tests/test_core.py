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
import numpy as np

class TestCore(unittest.TestCase):
    model_path = 'tests/fixtures/min_model.xlsx'
    target_specie_ids = ['RNA_1[c]','RNA_2[c]','RNA_3[c]','RNA_4[c]','RNA_5[c]']

    def test_is_constant_species(self):
        tweak_specie_ids = ['ATP[c]']
        is_constant, run_results = wc_test.ModelDynamicsTestCase().is_constant_species(model=self.model_path,
                  end_time=600, checkpoint_period=10, tweak_specie_ids=tweak_specie_ids, target_specie_ids=self.target_specie_ids)
        self.assertEqual(is_constant, ['False', 'False', 'False', 'True', 'True'])
        wc_lang.SpeciesType.objects.reset() # reset indexer

        tweak_specie_ids = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']
        is_constant, run_results = wc_test.ModelDynamicsTestCase().is_constant_species(model=self.model_path,
                  end_time=600, checkpoint_period=10, tweak_specie_ids=tweak_specie_ids, target_specie_ids=self.target_specie_ids)
        self.assertEqual(is_constant, ['False', 'False', 'False', 'True', 'True']) # Degradation active for first 3 species
        wc_lang.SpeciesType.objects.reset() # reset indexer

    def test_is_constant_reactions(self):
        tweak_reaction_ids = ['transcription_RNA_4']
        target_specie_ids = ['RNA_1[c]','RNA_2[c]','RNA_3[c]','RNA_4[c]','RNA_5[c]']
        is_constant, run_results = wc_test.ModelDynamicsTestCase().is_constant_reactions(model=self.model_path, end_time=600, checkpoint_period=10,
                  tweak_reaction_ids=tweak_reaction_ids, target_specie_ids=self.target_specie_ids)
        self.assertEqual(is_constant, ['False', 'False', 'False', 'True', 'True'])
        wc_lang.SpeciesType.objects.reset() # reset indexer

        tweak_reaction_ids = ['transcription_RNA_1', 'transcription_RNA_2', 'transcription_RNA_3', 'transcription_RNA_4', 'transcription_RNA_5']
        is_constant, run_results = wc_test.ModelDynamicsTestCase().is_constant_reactions(model=self.model_path, end_time=600, checkpoint_period=10,
                  tweak_reaction_ids=tweak_reaction_ids, target_specie_ids=self.target_specie_ids)
        self.assertEqual(is_constant, ['False', 'False', 'False', 'True', 'True']) # Degradation active for first 3 species
        wc_lang.SpeciesType.objects.reset() # reset indexer

    def test_scan_reactions(self):
        k_cat_vector = np.linspace(0,0.5,6)
        final_concentrations = wc_test.ModelDynamicsTestCase().scan_reactions(self.model_path, 900, 10, ['transcription_RNA_1'], 'RNA_1[c]', k_cat_vector)
        self.assertGreater(final_concentrations[1], final_concentrations[0])
        self.assertGreater(final_concentrations[2], final_concentrations[1])
        self.assertGreater(final_concentrations[3], final_concentrations[2])
        self.assertGreater(final_concentrations[4], final_concentrations[3])
        print('here')
        wc_lang.SpeciesType.objects.reset() # reset indexer

    def test_scan_species(self):
        init_concentrations = np.linspace(0,20000,6)
        final_concentrations = wc_test.ModelDynamicsTestCase().scan_species(self.model_path, 900, 10, ['RNA_4[c]'], 'RNA_4[c]', init_concentrations)
        self.assertGreater(final_concentrations[1], final_concentrations[0])
        self.assertGreater(final_concentrations[2], final_concentrations[1])
        self.assertGreater(final_concentrations[3], final_concentrations[2])
        self.assertGreater(final_concentrations[4], final_concentrations[3])
        wc_lang.SpeciesType.objects.reset() # reset indexer

    def test_reactions_mass_balanced(self):
        mass_balanced  = wc_test.ModelStaticTestCase().are_reactions_mass_balanced(self.model_path)
        # In test model transcription reactions are charge-, but not mass balanced;
        # degradation reactions are mass-, but not charge balanced
        self.assertEqual(mass_balanced['degradation_RNA_1'], True)
        self.assertEqual(mass_balanced['degradation_RNA_2'], True)
        self.assertEqual(mass_balanced['degradation_RNA_3'], True)
        self.assertEqual(mass_balanced['degradation_RNA_4'], True)
        self.assertEqual(mass_balanced['degradation_RNA_5'], True)
        self.assertEqual(mass_balanced['transcription_RNA_1'], False)
        self.assertEqual(mass_balanced['transcription_RNA_2'], False)
        self.assertEqual(mass_balanced['transcription_RNA_3'], False)
        self.assertEqual(mass_balanced['transcription_RNA_4'], False)
        self.assertEqual(mass_balanced['transcription_RNA_5'], False)
        wc_lang.SpeciesType.objects.reset()

    def test_reactions_charge_balanced(self):
        charge_balanced = wc_test.ModelStaticTestCase().are_reactions_charge_balanced(self.model_path)
        # In test model transcription reactions are charge-, but not mass balanced;
        # degradation reactions are mass-, but not charge balanced
        self.assertEqual(charge_balanced['degradation_RNA_1'], False)
        self.assertEqual(charge_balanced['degradation_RNA_2'], False)
        self.assertEqual(charge_balanced['degradation_RNA_3'], False)
        self.assertEqual(charge_balanced['degradation_RNA_4'], False)
        self.assertEqual(charge_balanced['degradation_RNA_5'], False)
        self.assertEqual(charge_balanced['transcription_RNA_1'], True)
        self.assertEqual(charge_balanced['transcription_RNA_2'], True)
        self.assertEqual(charge_balanced['transcription_RNA_3'], True)
        self.assertEqual(charge_balanced['transcription_RNA_4'], True)
        self.assertEqual(charge_balanced['transcription_RNA_5'], True)
