""" Methods for verifying models and simulations

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-05-10
:Copyright: 2018, Karr Lab
:License: MIT

TODO:
- all reaction methods: currently len(rate_laws)=1 assumed, generalize
- mod_parameters values are INTs in change_methods, but LISTs for sim_scan methods, synchornize
"""

import tempfile
import unittest
import wc_kb
import wc_kb.io
import wc_lang
import wc_lang.io
from wc_onto import onto
from wc_sim.simulation import Simulation
from wc_sim.run_results import RunResults


class KnowledgeBaseTestCase(unittest.TestCase):
    """ Methods for testing knowledge bases for WC models 

    Attributes:
        kb (:obj:`wc_kb.KnowledgeBase`): knowledge base

    Class attributes:
        KB (:obj:`wc_kb.KnowledgeBase` or :obj:`str`): knowledge base or path to a 
            knowledge base file
    """
    KB = None

    def setUp(self):
        if isinstance(self.KB, wc_kb.KnowledgeBase):
            self.kb = self.KB.copy()
        else:
            self.kb = wc_kb.io.Reader().run(self.KB)[wc_kb.KnowledgeBase][0]


class ModelTestCase(unittest.TestCase):
    """ Methods for testing WC models

    Attributes:        
        model (:obj:`wc_lang.Model`): model
        kb (:obj:`wc_kb.KnowledgeBase`): knowledge base
        results_dir (:obj:`str`): path to directory where results will be stored

    Class attributes:
        MODEL (:obj:`wc_lang.Model` or :obj:`str`): model or path to a model file
        KB (:obj:`wc_kb.KnowledgeBase` or :obj:`str`): knowledge base or path to a 
            knowledge base file
    """

    MODEL = None
    KB = None

    def setUp(self):
        if isinstance(self.MODEL, wc_lang.Model):
            self.model = self.MODEL.copy()
        else:
            self.model = wc_lang.io.Reader().run(self.MODEL)[wc_lang.Model][0]

        if isinstance(self.KB, wc_kb.KnowledgeBase):
            self.kb = self.KB.copy()
        elif self.KB is not None:
            self.kb = wc_kb.io.Reader().run(self.KB)[wc_kb.KnowledgeBase][0]

        self.results_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.results_dir)

    def get_species(self, id):
        return self.model.species.get_one(id=id)

    def get_reaction(self, id):
        return self.model.reactions.get_one(id=id)

    """ Methods to perturb model """
    # todo: use wc_lang.transform.ChangeValueTransform

    def select_submodels(self, mod_submodels):
        """ Turn off all submodels, except the ones listed in submodel_ids """
        for id, status in mod_submodels.items():
            if not status:
                submodel = self.model.submodels.get_one(id=id)
                for reaction in submodel.reactions:
                    for rate_law in reaction.rate_laws:
                        rate_law.expression.parameters.get_one(type=onto['WC:k_cat']).value = 0

    def change_parameter_values(self, mod_parameters):
        for id, value in mod_parameters.items():
            self.model.parameters.get_one(id=id).value = value

    def change_species_mean_init_concentrations(self, mod_species):
        for id, mean in mod_species.items():
            self.model.species.get_one(id=id).distribution_init_concentration.mean = mean

    def change_reaction_k_cat_parameter_values(self, mod_reactions):
        for id, k_cat_value in mod_reactions.items():
            reaction = self.model.reactions.get_one(id=id)
            reaction.rate_laws[0].expression.parameters.get_one(type=onto['WC:k_cat']).value = k_cat_value


class SimulationTestCase(ModelTestCase):
    """ Class to test simulations of models
    """

    """ Auxiliary methods """

    def simulate(self, end_time, checkpoint_period=None, n_sims=1):
        results = []

        simulation = Simulation(self.model)
        for i_sim in range(n_sims):
            temp_dir = tempfile.mkdtemp(dir=self.results_dir)
            num_events, results_dir = simulation.run(end_time=end_time,
                                                     results_dir=temp_dir,
                                                     checkpoint_period=checkpoint_period)

            run_results = RunResults(results_dir)
            results.append(run_results)
        return results

    """ Methods to obtain numbers to compare to exp data """

    def delta_conc(self, species, run_results):

        delta_conc = {}
        for specie_id in species:
            specie = self.get_species(specie_id)
            concentration = run_results.get('populations')[specie.id].values
            delta_conc[specie_id] = concentration[len(concentration)-1]-concentration[0]

        return delta_conc

    def avg_conc_time(self, target_specie_ids, end_time, checkpoint_period):
        # TODO: Test
        avg_conc = {}

        # Run model
        run_results = self.simulate(end_time=end_time, checkpoint_period=checkpoint_period)

        # Calculate avg concentration of target species
        for target_specie_id in target_specie_ids:
            target_specie = self.get_species(target_specie_id)
            concentration = run_results.get('populations')[target_specie.id][:].values  # convert panda.series to np.ndarray
            avg_conc[target_specie_id] = concentration.mean()

        return avg_conc

    def avg_conc_runs(self, n, target_specie_ids, end_time):
        # TODO: implement
        pass

    def get_growth_rate(self, end_time):
        # TODO: implement
        pass

    def sim_scan_parameters(self, mod_parameters, end_time, checkpoint_period):

        # Check if all dictionary values have same length
        lengths = []
        for parameter_id in mod_parameters:
            lengths.append(len(mod_parameters[parameter_id]))

        if lengths.count(lengths[0]) != len(lengths):
            raise SyntaxError('All values of mod_parameters should be a list with equal length')

        # Step over the defined parameter values and
        scan_results = []
        for value in mod_parameters[parameter_id]:  # step over each paramater value
            for parameter_id in mod_parameters:  # step over each parameters
                self.model.parameters.get_one(id=parameter_id).value = value

            run_result = self.simulate(end_time=end_time, checkpoint_period=checkpoint_period)[0]
            scan_results.append(run_result)

        return scan_results

    def sim_scan_species(self, mod_species, end_time, checkpoint_period):

        # Check if all dictionary values have same length
        lengths = []
        for species_id in mod_species:
            lengths.append(len(mod_species[species_id]))

        if lengths.count(lengths[0]) != len(lengths):
            raise SyntaxError('All values of mod_species should be a list with equal length')

        # Step over the defined parameter values and
        scan_results = []
        for value_index in range(0, lengths[0]):  # step over each paramater value
            for specie_id in mod_species:  # step over each parameters
                self.get_species(specie_id).distribution_init_concentration.mean = mod_species[specie_id][value_index]

            run_result = self.simulate(end_time=end_time, checkpoint_period=checkpoint_period)[0]
            scan_results.append(run_result)

        return scan_results

    def sim_scan_reactions(self, mod_reactions, end_time, checkpoint_period):
        # Check if all dictionary values have same length
        lengths = []
        for reactions_id in mod_reactions:
            lengths.append(len(mod_reactions[reactions_id]))

        if lengths.count(lengths[0]) != len(lengths):
            raise SyntaxError('All values of mod_reactions should be a list with equal length')

        # Step over the defined parameter values and
        scan_results = []
        for value_index in range(0, lengths[0]):  # step over each paramater value
            for reaction_id in mod_reactions:  # step over each parameters
                reaction = self.get_reaction(reaction_id)
                reaction.rate_laws[0].expression.parameters.get_one(
                    type=onto['WC:k_cat']).value = mod_reactions[reaction_id][value_index]

            run_result = self.simulate(end_time=end_time, checkpoint_period=checkpoint_period)[0]
            scan_results.append(run_result)

        return scan_results
