""" Test of wc_test.core

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-05-10
:Copyright: 2018, Karr Lab
:License: MIT

TODO:
- investigate what goes wrong with test_delta_conc when multiple species are passed as arguments
- checek whether multiple arguments in lists/dicts are testsed everywhere

"""

import wc_lang
import wc_test
import wc_sim
import unittest
import os
import numpy as np

class ModelTestCaseTests(unittest.TestCase):
    model_path = 'tests/fixtures/min_model.xlsx'
    model = wc_lang.io.Reader().run(model_path)

    def test_init(self):
        test_case1 = wc_test.ModelTestCase(model=self.model_path)
        test_case2 = wc_test.ModelTestCase(model=self.model, checkpoint_period=20, _results_dir='/home/test')

        self.assertIsInstance(test_case1.model, wc_lang.core.Model)
        self.assertEqual(test_case1.checkpoint_period, 30)
        self.assertEqual(test_case1._results_dir, os.path.expanduser('~/tmp/checkpoints_dir/'))

        self.assertIsInstance(test_case2.model, wc_lang.core.Model)
        self.assertEqual(test_case2.checkpoint_period, 20)
        self.assertEqual(test_case2._results_dir, '/home/test')

    def test_get_specie(self):
        specie = wc_test.ModelTestCase(model=self.model_path).get_specie('RNA_1[c]')
        self.assertIsInstance(specie, wc_lang.core.Species)
        self.assertEqual(specie.id(), 'RNA_1[c]')

        specie = wc_test.ModelTestCase(model=self.model_path).get_specie('C[c]')
        self.assertIsInstance(specie, wc_lang.core.Species)
        self.assertEqual(specie.id(), 'C[c]')

    def test_get_reaction(self):
        reaction = wc_test.ModelTestCase(model=self.model_path).get_reaction('transcription_RNA_1')
        self.assertIsInstance(reaction, wc_lang.core.Reaction)
        self.assertEqual(reaction.id, 'transcription_RNA_1')

        reaction = wc_test.ModelTestCase(model=self.model_path).get_reaction('degradation_RNA_1')
        self.assertIsInstance(reaction, wc_lang.core.Reaction)
        self.assertEqual(reaction.id, 'degradation_RNA_1')

class StaticTestCaseTests(unittest.TestCase):
    model_path = 'tests/fixtures/min_model.xlsx'
    # In test model transcription reactions are charge-, but not mass balanced;
    # degradation reactions are mass-, but not charge balanced

    def test_is_mass_balanced(self):
        self.assertTrue(wc_test.StaticTestCase(model=self.model_path).is_mass_balanced('degradation_RNA_1'))
        self.assertTrue(wc_test.StaticTestCase(model=self.model_path).is_mass_balanced('degradation_RNA_2'))
        self.assertFalse(wc_test.StaticTestCase(model=self.model_path).is_mass_balanced('transcription_RNA_1'))
        self.assertFalse(wc_test.StaticTestCase(model=self.model_path).is_mass_balanced('transcription_RNA_2'))

    def test_is_charge_balanced(self):
        self.assertFalse(wc_test.StaticTestCase(model=self.model_path).is_charge_balanced('degradation_RNA_1'))
        self.assertFalse(wc_test.StaticTestCase(model=self.model_path).is_charge_balanced('degradation_RNA_2'))
        self.assertTrue(wc_test.StaticTestCase(model=self.model_path).is_charge_balanced('transcription_RNA_1'))
        self.assertTrue(wc_test.StaticTestCase(model=self.model_path).is_charge_balanced('transcription_RNA_2'))

    def test_reactions_mass_balanced(self):
        mass_balanced  = wc_test.StaticTestCase(model=self.model_path).reactions_mass_balanced()
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
        charge_balanced = wc_test.StaticTestCase(self.model_path).reactions_charge_balanced()
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

class DynamicTestCaseTests(unittest.TestCase):
    model_path = 'tests/fixtures/min_model.xlsx'
    target_specie_ids = ['RNA_1[c]','RNA_2[c]','RNA_3[c]','RNA_4[c]','RNA_5[c]']

    def test_simulate(self):
        results = wc_test.DynamicTestCase(model=self.model_path).simulate(end_time=300)
        self.assertEqual(len(results), 1)
        self.assertIsInstance(results, list)
        self.assertIsInstance(results[0], wc_sim.multialgorithm.run_results.RunResults)

        results = wc_test.DynamicTestCase(model=self.model_path).simulate(end_time=2100, n=2)
        self.assertEqual(len(results), 2)
        self.assertIsInstance(results, list)
        self.assertIsInstance(results[0], wc_sim.multialgorithm.run_results.RunResults)
        self.assertIsInstance(results[1], wc_sim.multialgorithm.run_results.RunResults)

    def test_perturb_parameters(self):

        test_case = wc_test.DynamicTestCase(model=self.model_path)
        self.assertEqual(test_case.model.parameters.get_one(id='cell_cycle_length').value, 28800)
        self.assertEqual(test_case.model.parameters.get_one(id='fractionDryWeight').value, 0.7)

        mod_parameters={'cell_cycle_length':5,'fractionDryWeight':5}
        test_case.perturb_parameters(mod_parameters=mod_parameters)

        self.assertEqual(test_case.model.parameters.get_one(id='cell_cycle_length').value, 5)
        self.assertEqual(test_case.model.parameters.get_one(id='fractionDryWeight').value, 5)

    def test_perturb_species(self):

        test_case = wc_test.DynamicTestCase(model=self.model_path)
        self.assertEqual(test_case.get_specie('RNA_1[c]').concentration.value, 1000)
        self.assertEqual(test_case.get_specie('RNA_2[c]').concentration.value, 1000)

        mod_species={'RNA_1[c]':444,'RNA_2[c]':555}
        test_case.perturb_species(mod_species=mod_species)

        self.assertEqual(test_case.get_specie('RNA_1[c]').concentration.value, 444)
        self.assertEqual(test_case.get_specie('RNA_2[c]').concentration.value, 555)

    def test_perturb_reactions(self):

        test_case = wc_test.DynamicTestCase(model=self.model_path)
        self.assertEqual(test_case.get_reaction('transcription_RNA_1').rate_laws[0].k_cat, 0.05)
        self.assertEqual(test_case.get_reaction('degradation_RNA_1').rate_laws[0].k_cat, 0.035)

        mod_reactions={'transcription_RNA_1':5, 'degradation_RNA_1':6}
        test_case.perturb_reactions(mod_reactions=mod_reactions)

        self.assertEqual(test_case.get_reaction('transcription_RNA_1').rate_laws[0].k_cat, 5)
        self.assertEqual(test_case.get_reaction('degradation_RNA_1').rate_laws[0].k_cat, 6)

    def test_scan_parameters(self):
        """ Function to be built """
        pass

    def test_scan_species(self):
        """ Function to be built """
        pass

    def test_scan_reactions(self):
        """ Function to be built """
        pass

    def test_delta_conc(self):
        #test_case = wc_test.DynamicTestCase(model=self.model_path)
        #results = test_case.simulate(end_time=3000)
        #delta = test_case.delta_conc(['RNA_1[c],RNA_2[c]'], results[0])
        pass

    def test_avg_conc_time(self):
        """ Function to be built """
        pass

    def test_avg_conc_runs(self):
        """ Function to be built """
        pass

    def test_get_growth_rate(self):
        """ Function to be built """
        pass

# wc_lang.SpeciesType.objects.reset() # reset indexer
