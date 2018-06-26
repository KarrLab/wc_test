""" Methods for verifying models

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-05-10
:Copyright: 2018, Karr Lab
:License: MIT

TODO:
- generalize inputs; concentration / rate laws should be set to 0, but rather to a user defined value / list of values
- generalize test condition; concentration of target specie should not neccessarily stay contant, but rather stay within epsilon
"""

import os
import re
import wc_lang
import unittest
from wc_sim.multialgorithm.simulation import Simulation
from wc_sim.multialgorithm.run_results import RunResults

class KnowledgeBaseTestCase(unittest.TestCase):
    """ Test case for a knowledge base

    Attributes:
        kb (:obj:`wc_kb.KnowledgeBase`): knowledge base
    """
    pass

class ModelStaticTestCase(unittest.TestCase):
    """ Test case for static properties of models
    """

    def reactions_mass_balanced(self, model):
        """ Testing whether reactions in the model are mass balanced.

        Args:
            model (:obj:`wc_lang.core.Model`): model

        Returns:
            mass_balanced (:obj:`dict`): keys are reaction.ids; value is True if reaction is mass balanced, False otherwise
        """

        if not isinstance(model, wc_lang.core.Model):
            model = wc_lang.io.Reader().run(model)

        mass_balanced = {}
        imbalances = []

        for reaction in model.get_reactions():
            is_balanced = self.is_mass_balanced(reaction)
            mass_balanced[reaction.id] = is_balanced

            if not is_balanced:
                imbalances.append(reaction.id)

        if imbalances:
            print('The following reactions are not mass balanced:\n  {}'.format('\n  '.join(imbalances)))

        return mass_balanced

    def reactions_charge_balanced(self, model):
        """ Testing whether reactions in the model are charge balanced.

        Args:
            model (:obj:`wc_lang.core.Model`): model

        Returns:
            charge_balanced (:obj:`dict`): keys are reaction.ids; value is True if reaction is charge balanced, False otherwise
        """

        if not isinstance(model, wc_lang.core.Model):
            model = wc_lang.io.Reader().run(model)

        charge_balanced = {}
        imbalances = []

        for reaction in model.get_reactions():
            is_balanced = self.is_charge_balanced(reaction)
            charge_balanced[reaction.id] = is_balanced

            if not is_balanced:
                imbalances.append(reaction.id)

        if imbalances:
            print('The following reactions are not charge balanced:\n  {}'.format('\n  '.join(imbalances)))

        return charge_balanced

    def is_charge_balanced(self, reaction):
        """ Testing whether a reaction is charge balanced. Retruns boolean True if it is, False otherwise

        Args:
            reaction (:obj:`wc_lang.core.Reaction`): reaction
        """

        if not isinstance(reaction, wc_lang.core.Reaction):
            raise Exception('{} is not an instance of the expected wc_lang.core.reaction class'.format(reaction))

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

    def is_mass_balanced(self, reaction):
        """ Testing whether a reaction is mass balanced. Retruns boolean True if it is, False otherwise

        Args:
            reaction (:obj:`wc_lang.core.Reaction`): reaction
        """

        if not isinstance(reaction, wc_lang.core.Reaction):
            raise Exception('{} is not an instance of the expected wc_lang.core.reaction class'.format(reaction))

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

class ModelDynamicsTestCase(unittest.TestCase):
    """ Class to test dynamic properties of models

        Attributes:
            model (:obj:`wc_lang.core.Model`) OR (:obj:`str`): model or path to the model file
            checkpoint_period (:obj:`int`): interval at which results are saved
            results_dir (:obj:`str`): path to directory where results will be stored
    """

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

    def setUp(self):
        self._results_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._results_dir)

    def simulate(self, end_time, n=None): #NEED TO ADD TEST!

        if not n:
            n=1

        model = self.model
        checkpoint_period = self.checkpoint_period
        results_dir = self._results_dir

        simulation = Simulation(model)
        for n in range(0,n):
            num_events, results_dir = simulation.run(end_time=end_time, results_dir=results_dir, checkpoint_period=checkpoint_period)

        run_results = RunResults(results_dir)
        return run_results

    def get_specie(self, specie_id): #NEED TO ADD TEST!
        model = self.model

        temp = re.findall('[a-zA-Z_0-9]+', specie_id)
        specie_type_id = temp[0]
        compartment_id = temp[1]
        specie_compartment = model.compartments.get_one(id=compartment_id)
        specie = model.species_types.get_one(id=specie_type_id).species.get_one(compartment=specie_compartment)

        return specie

    def is_constant_species(self, tweak_specie_ids, target_specie_ids, end_time):
        """ Checks whether setting the concentration(s) of tweak_specie_ids to 0, stabilizes (i.e. remains constant) the concentration of
            target_specie_ids within end_time

            Args:
                species_type (:obj:`wc_lang.SpeciesType`): species type
                target_species_type (:obj:`wc_lang.SpeciesType`):
                time (:obj:`float`): time

            Returns:
                is_constant (:obj:`list`): list of boolean True and False
                run_results (:obj:`wc_sim.multialgorithm.run_results.RunResults`): run results
        """

        # Set concentration of species to 0:
        for tweak_specie_id in tweak_specie_ids:
            tweak_specie = self.get_specie(tweak_specie_id)
            tweak_specie.concentration.value = 0

        # Run model
        run_results = self.simulate(end_time=end_time)

        # Check target species concentration at the end of simulation
        is_constant = []
        for target_specie_id in target_specie_ids:
            target_specie = self.get_specie(target_specie_id)
            concentration = run_results.get('populations')[target_specie.id()][:].values #convert panda.series to np.ndarray

            if all(concentration[0]==concentration):
                is_constant.append('True')
            else:
                is_constant.append('False')

        return is_constant

    def scan_species(self, tweak_specie_ids, target_specie_id, init_concentrations, end_time):

        final_concentrations=[]
        target_specie = self.get_specie(target_specie_id)

        for init_concentration in init_concentrations:
            for tweak_specie_id in tweak_specie_ids:
                tweak_specie = self.get_specie(tweak_specie_id)
                tweak_specie.concentration.value = init_concentration

            # Run model
            run_results = self.simulate(end_time=end_time)

            concentrations = run_results.get('populations')[target_specie.id()][:].values
            final_time = int(end_time/self.checkpoint_period)
            concentration = concentrations[final_time]
            final_concentrations.append(concentration)

        return final_concentrations

    def is_constant_reactions(self, tweak_reaction_ids, target_specie_ids, end_time):

        # Set concentration of species to 0:
        for tweak_reaction_id in tweak_reaction_ids:
            for reaction in self.model.get_reactions():
                if reaction.id == tweak_reaction_id:
                    reaction.rate_laws[0].k_m = 0
                    reaction.rate_laws[0].k_cat = 0
                    break

        # Run model
        run_results = self.simulate(end_time=end_time)

        # Check target species concentration at the end of simulation
        is_constant = []
        for target_specie_id in target_specie_ids:
            target_specie = self.get_specie(target_specie_id)
            concentration = run_results.get('populations')[target_specie.id()][:].values #convert panda.series to np.ndarray

            if all(concentration[0]==concentration):
                is_constant.append('True')
            else:
                is_constant.append('False')

        return is_constant

    def scan_reactions(self, tweak_reaction_ids, target_specie_id, k_cats, end_time ):
        final_concentrations=[]
        target_specie = self.get_specie(target_specie_id)

        for k_cat in k_cats:
            # Set concentration of species to 0:
            for tweak_reaction_id in tweak_reaction_ids:
                for reaction in self.model.get_reactions():
                    if reaction.id == tweak_reaction_id:
                        reaction.rate_laws[0].k_cat = k_cat
                        break

            # Run model
            run_results = self.simulate(end_time=end_time)

            concentrations = run_results.get('populations')[target_specie.id()][:].values
            final_time = int(end_time/self.checkpoint_period)
            concentration = concentrations[final_time]
            final_concentrations.append(concentration)

        return final_concentrations

    def avg_conc_time(self, target_specie_ids, end_time): #NEED TO ADD TEST!
        avg_conc = {}

        # Run model
        run_results = self.simulate(end_time=end_time)

        # Calculate avg concentration of target species
        for target_specie_id in target_specie_ids:
            target_specie = self.get_specie(target_specie_id)
            concentration = run_results.get('populations')[target_specie.id()][:].values #convert panda.series to np.ndarray
            avg_conc[target_specie_id] = concentration.mean()

        return avg_conc

    def avg_conc_runs(self, n, target_specie_ids, end_time):
        pass

    def parameter_scan(self, parameter, end_time):
        pass

    def get_growth_rate(self, end_time):
        pass
