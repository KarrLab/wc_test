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

    def get_specie(self, specie_id):
        temp = re.findall('[a-zA-Z_0-9]+', specie_id)
        specie_type_id = temp[0]
        compartment_id = temp[1]
        specie_compartment = self.model.compartments.get_one(id=compartment_id)
        specie = self.model.species_types.get_one(id=specie_type_id).species.get_one(compartment=specie_compartment)

        return specie

    def get_reaction(self, reaction_id):
        for reaction in self.model.get_reactions():
            if reaction.id == reaction_id:
                return reaction

    def select_submodels(self, mod_submodels):
        """ Turn off all submodels, except the ones listed in submodel_ids """

        for submodel in self.model.submodels:
            if mod_submodels[submodel.id] == True:
                continue
            else:
                for reaction in submodel.reactions:
                    reaction.rate_laws[0].k_cat=0

    """ Methods to perturb model """
    def perturb_parameters(self, mod_parameters):

        for parameter_id in mod_parameters:
            self.model.parameters.get_one(id=parameter_id).value = mod_parameters[parameter_id]

    def perturb_species(self, mod_species):

        for specie_id in mod_species:
            self.get_specie(specie_id).concentration.value = mod_species[specie_id]

    def perturb_reactions(self, mod_reactions):
        for reaction_id in mod_reactions:
            reaction = self.get_reaction(reaction_id)
            reaction.rate_laws[0].k_cat = mod_reactions[reaction_id]

class StaticTestCase(ModelTestCase):
    """ Test case for static properties of models """

    def check_init_compartment_volumes(self, bounds):
        check_bounds={}
        for compartment in self.model.compartments:
            if bounds[0] < compartment.initial_volume < bounds[1]:
                check_bounds[compartment.id]=True
            else:
                check_bounds[compartment.id]=False
        return check_bounds

    def check_init_species_types_charges(self, bounds):
        check_bounds={}
        for species_type in self.model.species_types:
            if bounds[0] <= species_type.charge <= bounds[1]:
                check_bounds[species_type.id]=True
            else:
                check_bounds[species_type.id]=False
        return check_bounds

    def check_init_species_types_weights(self, bounds):
        check_bounds={}
        for species_type in self.model.species_types:
            if bounds[0] <= species_type.molecular_weight <= bounds[1]:
                check_bounds[species_type.id] = True
            else:
                check_bounds[species_type.id] = False
        return check_bounds

    def check_init_reactions_rates(self, bounds):
        check_bounds={}
        for reaction in self.model.get_reactions():
            if bounds[0] <= reaction.rate_laws[0].k_cat <= bounds[1]:
                check_bounds[reaction.id] = True
            else:
                check_bounds[reaction.id] = False
        return check_bounds

    def reactions_mass_balanced(self):
        """ Testing whether reactions in the model are mass balanced.

        Args:
            model (:obj:`wc_lang.core.Model`): model

        Returns:
            mass_balanced (:obj:`dict`): keys are reaction.ids; value is True if reaction is mass balanced, False otherwise
        """

        mass_balanced = {}
        imbalances = []

        for reaction in self.model.get_reactions():
            is_balanced = self.is_mass_balanced(reaction.id)
            mass_balanced[reaction.id] = is_balanced

            if not is_balanced:
                imbalances.append(reaction.id)

        if imbalances:
            print('The following reactions are not mass balanced:\n  {}'.format('\n  '.join(imbalances)))

        return mass_balanced

    def reactions_charge_balanced(self):
        """ Testing whether reactions in the model are charge balanced.

        Args:
            model (:obj:`wc_lang.core.Model`): model

        Returns:
            charge_balanced (:obj:`dict`): keys are reaction.ids; value is True if reaction is charge balanced, False otherwise
        """

        charge_balanced = {}
        imbalances = []

        for reaction in self.model.get_reactions():
            is_balanced = self.is_charge_balanced(reaction.id)
            charge_balanced[reaction.id] = is_balanced

            if not is_balanced:
                imbalances.append(reaction.id)

        if imbalances:
            print('The following reactions are not charge balanced:\n  {}'.format('\n  '.join(imbalances)))

        return charge_balanced

    def is_charge_balanced(self, reaction_id):
        """ Testing whether a reaction is charge balanced. Retruns boolean True if it is, False otherwise

        Args:
            reaction (:obj:`wc_lang.core.Reaction`): reaction
        """

        reaction = self.get_reaction(reaction_id)
        lhs_charge = 0
        rhs_charge = 0

        for participant in reaction.participants:
            if participant.coefficient < 0:
                lhs_charge += abs(participant.coefficient)*participant.species.species_type.charge
            elif participant.coefficient > 0:
                rhs_charge += abs(participant.coefficient)*participant.species.species_type.charge

        if lhs_charge == rhs_charge:
            return True
        else:
            return False

    def is_mass_balanced(self, reaction_id):
        """ Testing whether a reaction is mass balanced. Retruns boolean True if it is, False otherwise

        Args:
            reaction (:obj:`wc_lang.core.Reaction`): reaction
        """

        reaction = self.get_reaction(reaction_id)
        lhs_mass = 0
        rhs_mass = 0

        for participant in reaction.participants:
            if participant.coefficient < 0:
                lhs_mass += abs(participant.coefficient)*participant.species.species_type.molecular_weight
            elif participant.coefficient > 0:
                rhs_mass += abs(participant.coefficient)*participant.species.species_type.molecular_weight

        if lhs_mass == rhs_mass:
            return True
        else:
            return False

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

        results=[]
        if not n:
            n=1

        simulation = Simulation(self.model)
        for n in range(0,n):
            num_events, results_dir = simulation.run(end_time = end_time,
                                                     results_dir = self._results_dir,
                                                     checkpoint_period = self.checkpoint_period)

            run_results = RunResults(results_dir)
            results.append(run_results)
        return results

    """ Methods to obtain numbers to compare to exp data """
    def delta_conc(self, species, run_results):

        delta_conc = {}
        for specie_id in species:
            specie = self.get_specie(specie_id)
            concentration = run_results.get('populations')[specie.id].values
            delta_conc[specie_id] = concentration[len(concentration)-1]-concentration[0]

        return delta_conc

    def avg_conc_time(self, target_specie_ids, end_time): #NEED TO ADD TEST!
        avg_conc = {}

        # Run model
        run_results = self.simulate(end_time=end_time)

        # Calculate avg concentration of target species
        for target_specie_id in target_specie_ids:
            target_specie = self.get_specie(target_specie_id)
            concentration = run_results.get('populations')[target_specie.id][:].values #convert panda.series to np.ndarray
            avg_conc[target_specie_id] = concentration.mean()

        return avg_conc

    def avg_conc_runs(self, n, target_specie_ids, end_time):
        pass

    def get_growth_rate(self, end_time):
        pass

    def sim_scan_parameters(self, mod_parameters, end_time):

        # Check if all dictionary values have same length
        lengths=[]
        for parameter_id in mod_parameters:
             lengths.append(len(mod_parameters[parameter_id]))

        if lengths.count(lengths[0]) != len(lengths):
            raise SyntaxError('All values of mod_parameters should be a list with equal length')

        # Step over the defined parameter values and
        scan_results = []
        for value_index in range(0, lengths[0]): # step over each paramater value
            for parameter_id in mod_parameters: # step over each parameters
                self.model.parameters.get_one(id=parameter_id).value = mod_parameters[parameter_id][value_index]

            run_result = self.simulate(end_time=end_time)[0]
            scan_results.append(run_result)

        return scan_results

    def sim_scan_species(self, mod_species, end_time):

        # Check if all dictionary values have same length
        lengths=[]
        for species_id in mod_species:
             lengths.append(len(mod_species[species_id]))

        if lengths.count(lengths[0]) != len(lengths):
            raise SyntaxError('All values of mod_species should be a list with equal length')

        # Step over the defined parameter values and
        scan_results = []
        for value_index in range(0, lengths[0]): # step over each paramater value
            for specie_id in mod_species: # step over each parameters
                self.get_specie(specie_id).concentration.value = mod_species[specie_id][value_index]

            run_result = self.simulate(end_time=end_time)[0]
            scan_results.append(run_result)

        return scan_results

    def sim_scan_reactions(self, mod_reactions, end_time):
        # Check if all dictionary values have same length
        lengths=[]
        for reactions_id in mod_reactions:
             lengths.append(len(mod_reactions[reactions_id]))

        if lengths.count(lengths[0]) != len(lengths):
            raise SyntaxError('All values of mod_reactions should be a list with equal length')

        # Step over the defined parameter values and
        scan_results = []
        for value_index in range(0, lengths[0]): # step over each paramater value
            for reaction_id in mod_reactions: # step over each parameters
                reaction = self.get_reaction(reaction_id)
                reaction.rate_laws[0].k_cat = mod_reactions[reaction_id][value_index]

            run_result = self.simulate(end_time=end_time)[0]
            scan_results.append(run_result)

        return scan_results
