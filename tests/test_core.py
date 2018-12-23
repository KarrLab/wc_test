""" Test of wc_test.core

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-05-10
:Copyright: 2018, Karr Lab
:License: MIT

TODO:
- investigate what goes wrong with test_delta_conc when multiple species are passed as arguments
- checek whether multiple arguments in lists/dicts are testsed everywhere

"""

import os
import unittest
import wc_kb.io
import wc_lang
import wc_lang.io
import wc_sim
import wc_test.core


class KnowledgeBaseTestCaseTestCase(unittest.TestCase):
    KB_PATH = 'tests/fixtures/min_kb.xlsx'

    @unittest.skip('Todo')
    def setUp(self):
        self.kb = wc_kb.io.Reader().run(self.KB_PATH)

        class TestCase(wc_test.core.KnowledgeBaseTestCase):
            KB = self.kb

        self.test_case = TestCase()
        self.test_case.setUp()


class ModelTestCaseTestCase(unittest.TestCase):
    MODEL_PATH = 'tests/fixtures/min_model.xlsx'

    def setUp(self):
        self.model = wc_lang.io.Reader().run(self.MODEL_PATH)

        class TestCase(wc_test.core.ModelTestCase):
            MODEL = self.model

        self.test_case = TestCase()
        self.test_case.setUp()

    def test_setUp(self):
        self.assertTrue(self.test_case.model.is_equal(self.model))
        self.assertTrue(os.path.isdir(self.test_case.results_dir))

    def test_select_submodels(self):
        test_case = self.test_case
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
        species = self.test_case.get_species('RNA_1[c]')
        self.assertIsInstance(species, wc_lang.core.Species)
        self.assertEqual(species.id, 'RNA_1[c]')

        species = self.test_case.get_species('H[c]')
        self.assertIsInstance(species, wc_lang.core.Species)
        self.assertEqual(species.id, 'H[c]')

    def test_get_reaction(self):
        reaction = self.test_case.get_reaction('transcription_RNA_1')
        self.assertIsInstance(reaction, wc_lang.core.Reaction)
        self.assertEqual(reaction.id, 'transcription_RNA_1')

        reaction = self.test_case.get_reaction('degradation_RNA_1')
        self.assertIsInstance(reaction, wc_lang.core.Reaction)
        self.assertEqual(reaction.id, 'degradation_RNA_1')

    def test_change_parameter_values(self):
        test_case = self.test_case
        self.assertEqual(test_case.model.parameters.get_one(id='mean_doubling_time').value, 28800)

        mod_parameters = {'mean_doubling_time': 5}
        test_case.change_parameter_values(mod_parameters=mod_parameters)

        self.assertEqual(test_case.model.parameters.get_one(id='mean_doubling_time').value, 5)

    def test_change_species_mean_init_concentrations(self):
        test_case = self.test_case
        self.assertEqual(test_case.get_species('RNA_1[c]').distribution_init_concentration.mean, 1000)
        self.assertEqual(test_case.get_species('RNA_2[c]').distribution_init_concentration.mean, 1000)

        mod_species = {'RNA_1[c]': 444, 'RNA_2[c]': 555}
        test_case.change_species_mean_init_concentrations(mod_species=mod_species)

        self.assertEqual(test_case.get_species('RNA_1[c]').distribution_init_concentration.mean, 444)
        self.assertEqual(test_case.get_species('RNA_2[c]').distribution_init_concentration.mean, 555)

    def test_change_reaction_k_cat_parameter_values(self):
        test_case = self.test_case
        self.assertEqual(test_case.get_reaction('transcription_RNA_1').rate_laws[0].
                         expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value, 0.05)
        self.assertEqual(test_case.get_reaction('degradation_RNA_1').rate_laws[0].
                         expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value, 0.035)

        mod_reactions = {'transcription_RNA_1': 5, 'degradation_RNA_1': 6}
        test_case.change_reaction_k_cat_parameter_values(mod_reactions=mod_reactions)

        self.assertEqual(test_case.get_reaction('transcription_RNA_1').rate_laws[0].
                         expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value, 5)
        self.assertEqual(test_case.get_reaction('degradation_RNA_1').rate_laws[0].
                         expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value, 6)


class SimulationTestCaseTestCase(unittest.TestCase):
    MODEL_PATH = 'tests/fixtures/min_model.xlsx'

    def setUp(self):
        self.model = wc_lang.io.Reader().run(self.MODEL_PATH)

        class TestCase(wc_test.core.SimulationTestCase):
            MODEL = self.model

        self.test_case = TestCase()
        self.test_case.setUp()

    def test_simulate(self):
        results = self.test_case.simulate(end_time=10., checkpoint_period=5.)
        self.assertEqual(len(results), 1)
        self.assertIsInstance(results, list)
        self.assertIsInstance(results[0], wc_sim.multialgorithm.run_results.RunResults)

        results = self.test_case.simulate(end_time=10., checkpoint_period=5., n_sims=2)
        self.assertEqual(len(results), 2)
        self.assertIsInstance(results, list)
        self.assertIsInstance(results[0], wc_sim.multialgorithm.run_results.RunResults)
        self.assertIsInstance(results[1], wc_sim.multialgorithm.run_results.RunResults)

    @unittest.skip('Todo: implement')
    def test_delta_conc(self):
        #test_case = self.test_case
        #results = test_case.simulate(end_time=10., checkpoint_period=5.)
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
        test_case = self.test_case

        mod_parameters = {'mean_doubling_time': [5]}
        scan_results = test_case.sim_scan_parameters(mod_parameters=mod_parameters, end_time=10., checkpoint_period=5.)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

        mod_parameters = {'mean_doubling_time': [5, 6, 7]}
        scan_results = test_case.sim_scan_parameters(mod_parameters=mod_parameters, end_time=10., checkpoint_period=5.)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

    def test_sim_scan_species(self):
        test_case = self.test_case

        mod_species = {'RNA_1[c]': [5]}
        scan_results = test_case.sim_scan_species(mod_species=mod_species, end_time=10., checkpoint_period=5.)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

        mod_parameters = {'RNA_1[c]': [5, 6, 7], 'RNA_2[c]': [8, 9, 10]}
        scan_results = test_case.sim_scan_species(mod_species=mod_species, end_time=10., checkpoint_period=5.)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

    def test_sim_scan_reactions(self):
        test_case = self.test_case

        mod_reactions = {'transcription_RNA_1': [0.65]}
        scan_results = test_case.sim_scan_reactions(mod_reactions=mod_reactions, end_time=10., checkpoint_period=5.)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)

        mod_reactions = {'transcription_RNA_1': [0.55, 0.65, 0.75], 'degradation_RNA_1': [0.35, 0.45, 0.55]}
        scan_results = test_case.sim_scan_reactions(mod_reactions=mod_reactions, end_time=10., checkpoint_period=5.)
        self.assertIsInstance(scan_results, list)
        self.assertIsInstance(scan_results[0], wc_sim.multialgorithm.run_results.RunResults)
