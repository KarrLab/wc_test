""" Methods for verifying models

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-05-10
:Copyright: 2018, Karr Lab
:License: MIT

TODO:
- all reaction methods: currently len(rate_laws)=1 assumed, generalize
- mod_parameters values are INTs in perturb_methods, but LISTs for sim_scan methods, synchornize
"""

import os
import re
import wc_lang
import unittest
from wc_sim.multialgorithm.simulation import Simulation
from wc_sim.multialgorithm.run_results import RunResults


class ModelTestCase(unittest.TestCase):
    """ Base classe for WC model testing classes """

    def __init__(self, model, checkpoint_period=None, _results_dir=None):
        if not isinstance(model, wc_lang.core.Model):
            model = wc_lang.io.Reader().run(model)

        if not _results_dir:
            _results_dir = os.path.expanduser('~/tmp/checkpoints_dir/')

        if not checkpoint_period:
            checkpoint_period = 30

        self.model = model
        self._results_dir = _results_dir
        self.checkpoint_period = checkpoint_period

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
                        rate_law.expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value = 0

    def perturb_parameter_values(self, mod_parameters):
        for id, value in mod_parameters.items():
            self.model.parameters.get_one(id=id).value = value

    def perturb_species_mean_init_concentrations(self, mod_species):
        for id, mean in mod_species.items():
            self.model.species.get_one(id=id).distribution_init_concentration.mean = mean

    def perturb_reaction_k_cat_parameter_values(self, mod_reactions):
        for id, k_cat_value in mod_reactions.items():
            reaction = self.model.reactions.get_one(id=id)
            reaction.rate_laws[0].expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value = k_cat_value


class StaticTestCase(ModelTestCase):
    """ Test case for static properties of models """
    # TODO: implement meaningful tests
    pass


class DynamicTestCase(ModelTestCase):
    """ Class to test dynamic properties of models

        Attributes:
            model (:obj:`wc_lang.core.Model`) OR (:obj:`str`): model or path to the model file
            checkpoint_period (:obj:`int`): interval at which results are saved
            results_dir (:obj:`str`): path to directory where results will be stored
    """

    """ Auxiliary methods """

    def setUp(self):
        self._results_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._results_dir)

    def simulate(self, end_time, n=None):

        results = []
        if not n:
            n = 1

        simulation = Simulation(self.model)
        for n in range(0, n):
            num_events, results_dir = simulation.run(end_time=end_time,
                                                     results_dir=self._results_dir,
                                                     checkpoint_period=self.checkpoint_period)

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

    def avg_conc_time(self, target_specie_ids, end_time):
        # TODO: Test
        avg_conc = {}

        # Run model
        run_results = self.simulate(end_time=end_time)

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

    def sim_scan_parameters(self, mod_parameters, end_time):

        # Check if all dictionary values have same length
        lengths = []
        for parameter_id in mod_parameters:
            lengths.append(len(mod_parameters[parameter_id]))

        if lengths.count(lengths[0]) != len(lengths):
            raise SyntaxError('All values of mod_parameters should be a list with equal length')

        # Step over the defined parameter values and
        scan_results = []
        for value_index in range(0, lengths[0]):  # step over each paramater value
            for parameter_id in mod_parameters:  # step over each parameters
                self.model.parameters.get_one(id=parameter_id).value = mod_parameters[parameter_id][value_index]

            run_result = self.simulate(end_time=end_time)[0]
            scan_results.append(run_result)

        return scan_results

    def sim_scan_species(self, mod_species, end_time):

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

            run_result = self.simulate(end_time=end_time)[0]
            scan_results.append(run_result)

        return scan_results

    def sim_scan_reactions(self, mod_reactions, end_time):
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
                    type=wc_lang.ParameterType.k_cat).value = mod_reactions[reaction_id][value_index]

            run_result = self.simulate(end_time=end_time)[0]
            scan_results.append(run_result)

        return scan_results
