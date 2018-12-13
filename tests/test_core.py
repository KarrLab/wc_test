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
    MODEL_PATH = 'tests/fixtures/min_model.xlsx'

    def setUp(self):
        self.model = wc_lang.io.Reader().run(self.MODEL_PATH)

    def test_init(self):
        test_case1 = wc_test.ModelTestCase(model=self.MODEL_PATH)
        test_case2 = wc_test.ModelTestCase(model=self.model, checkpoint_period=20, _results_dir='/home/test')

        self.assertIsInstance(test_case1.model, wc_lang.core.Model)
        self.assertEqual(test_case1.checkpoint_period, 30)
        self.assertEqual(test_case1._results_dir, os.path.expanduser('~/tmp/checkpoints_dir/'))

        self.assertIsInstance(test_case2.model, wc_lang.core.Model)
        self.assertEqual(test_case2.checkpoint_period, 20)
        self.assertEqual(test_case2._results_dir, '/home/test')

    def test_select_submodels(self):

        test_case = wc_test.DynamicTestCase(model=self.MODEL_PATH)
        mod_submodels = {
            'transcription': True,
            'degradation': False,
        }
        test_case.select_submodels(mod_submodels=mod_submodels)

        for reaction in test_case.model.get_reactions():
            if mod_submodels[reaction.submodel.id] == False:
                print(reaction.submodel.id)
                self.assertEqual(reaction.rate_laws[0].
                                 expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value, 0)

    def test_get_species(self):
        species = wc_test.ModelTestCase(model=self.MODEL_PATH).get_species('RNA_1[c]')
        self.assertIsInstance(species, wc_lang.core.Species)
        self.assertEqual(species.id, 'RNA_1[c]')

        species = wc_test.ModelTestCase(model=self.MODEL_PATH).get_species('H[c]')
        self.assertIsInstance(species, wc_lang.core.Species)
        self.assertEqual(species.id, 'H[c]')

    def test_get_reaction(self):
        reaction = wc_test.ModelTestCase(model=self.MODEL_PATH).get_reaction('transcription_RNA_1')
        self.assertIsInstance(reaction, wc_lang.core.Reaction)
        self.assertEqual(reaction.id, 'transcription_RNA_1')

        reaction = wc_test.ModelTestCase(model=self.MODEL_PATH).get_reaction('degradation_RNA_1')
        self.assertIsInstance(reaction, wc_lang.core.Reaction)
        self.assertEqual(reaction.id, 'degradation_RNA_1')

    def test_perturb_parameter_values(self):
        test_case = wc_test.ModelTestCase(model=self.MODEL_PATH)
        self.assertEqual(test_case.model.parameters.get_one(id='cell_cycle_length').value, 28800)
        self.assertEqual(test_case.model.parameters.get_one(id='fractionDryWeight').value, 0.7)

        mod_parameters = {'cell_cycle_length': 5, 'fractionDryWeight': 5}
        test_case.perturb_parameter_values(mod_parameters=mod_parameters)

        self.assertEqual(test_case.model.parameters.get_one(id='cell_cycle_length').value, 5)
        self.assertEqual(test_case.model.parameters.get_one(id='fractionDryWeight').value, 5)

    def test_perturb_species_mean_init_concentrations(self):
        test_case = wc_test.ModelTestCase(model=self.MODEL_PATH)
        self.assertEqual(test_case.get_species('RNA_1[c]').distribution_init_concentration.mean, 1000)
        self.assertEqual(test_case.get_species('RNA_2[c]').distribution_init_concentration.mean, 1000)

        mod_species = {'RNA_1[c]': 444, 'RNA_2[c]': 555}
        test_case.perturb_species_mean_init_concentrations(mod_species=mod_species)

        self.assertEqual(test_case.get_species('RNA_1[c]').distribution_init_concentration.mean, 444)
        self.assertEqual(test_case.get_species('RNA_2[c]').distribution_init_concentration.mean, 555)

    def test_perturb_reaction_k_cat_parameter_values(self):
        test_case = wc_test.ModelTestCase(model=self.MODEL_PATH)
        self.assertEqual(test_case.get_reaction('transcription_RNA_1').rate_laws[0].
                         expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value, 0.05)
        self.assertEqual(test_case.get_reaction('degradation_RNA_1').rate_laws[0].
                         expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value, 0.035)

        mod_reactions = {'transcription_RNA_1': 5, 'degradation_RNA_1': 6}
        test_case.perturb_reaction_k_cat_parameter_values(mod_reactions=mod_reactions)

        self.assertEqual(test_case.get_reaction('transcription_RNA_1').rate_laws[0].
                         expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value, 5)
        self.assertEqual(test_case.get_reaction('degradation_RNA_1').rate_laws[0].
                         expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value, 6)


class StaticTestCaseTests(ModelTestCaseTests):
    pass


class DynamicTestCaseTests(ModelTestCaseTests):
    MODEL_PATH = 'tests/fixtures/min_model.xlsx'

    def test_simulate(self):
        results = wc_test.DynamicTestCase(model=self.MODEL_PATH).simulate(end_time=300)
        self.assertEqual(len(results), 1)
        self.assertIsInstance(results, list)
        self.assertIsInstance(results[0], wc_sim.multialgorithm.run_results.RunResults)

        results = wc_test.DynamicTestCase(model=self.MODEL_PATH).simulate(end_time=2100, n=2)
        self.assertEqual(len(results), 2)
        self.assertIsInstance(results, list)
        self.assertIsInstance(results[0], wc_sim.multialgorithm.run_results.RunResults)
        self.assertIsInstance(results[1], wc_sim.multialgorithm.run_results.RunResults)

    @unittest.skip('Todo: implement')
    def test_delta_conc(self):
        #test_case = wc_test.DynamicTestCase(model=self.MODEL_PATH)
        #results = test_case.simulate(end_time=3000)
        #delta = test_case.delta_conc(['RNA_1[c],RNA_2[c]'], results[0])
        pass

    @unittest.skip('Todo: implement')
    def test_avg_conc_time(self):
        pass

    @unittest.skip('Todo: implement')
    def test_avg_conc_runs(self):
        pass

    @unittest.skip('Todo: implement')
    def test_get_growth_rate(self):
        pass

    def test_sim_scan_parameters(self):
        test_case = wc_test.DynamicTestCase(model=self.MODEL_PATH)

        mod_parameters = {'cell_cycle_length': [5], 'fractionDryWeight': [5]}
        scan_results = test_case.sim_scan_parameters(mod_parameters=mod_parameters, end_time=600)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

        mod_parameters = {'cell_cycle_length': [5, 6, 7], 'fractionDryWeight': [8, 9, 10]}
        scan_results = test_case.sim_scan_parameters(mod_parameters=mod_parameters, end_time=600)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

    def test_sim_scan_species(self):
        test_case = wc_test.DynamicTestCase(model=self.MODEL_PATH)

        mod_species = {'RNA_1[c]': [5]}
        scan_results = test_case.sim_scan_species(mod_species=mod_species, end_time=600)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

        mod_parameters = {'RNA_1[c]': [5, 6, 7], 'RNA_2[c]': [8, 9, 10]}
        scan_results = test_case.sim_scan_species(mod_species=mod_species, end_time=600)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

    def test_sim_scan_reactions(self):
        test_case = wc_test.DynamicTestCase(model=self.MODEL_PATH)

        mod_reactions = {'transcription_RNA_1': [0.65]}
        scan_results = test_case.sim_scan_reactions(mod_reactions=mod_reactions, end_time=600)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

        mod_reactions = {'transcription_RNA_1': [0.55, 0.65, 0.75], 'degradation_RNA_1': [0.35, 0.45, 0.55]}
        scan_results = test_case.sim_scan_reactions(mod_reactions=mod_reactions, end_time=600)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)
