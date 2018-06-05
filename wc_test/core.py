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

class SubmodelDynamicsTestCase(unittest.TestCase):
    """ Methods for verifying submodels  """

    def is_constant_tweak_species(self, model_path, end_time, checkpoint_period, tweak_specie_ids, target_specie_ids):
        """ Checks whether setting the concentration of species_type[compartment] to 0 either:

            * reduces the concentration of specie 'target_id' (instance of 'target_class') to 0 within 'time'
            * stabilizes (i.e. remains constant) the concentration of target_species_type[target_compartment] within 'time'

            Args:
                model (:obj:`wc_lang.Model`): model
                species_type (:obj:`wc_lang.SpeciesType`): species type
                compartment (:obj:`wc_lang.CompartmentModel`): compartment
                target_species_type (:obj:`wc_lang.SpeciesType`):
                target_compartment (:obj:`wc_lang.CompartmentModel`):
                time (:obj:`float`): time
                is_constant (:obj:`bool`):

            Returns:
                :obj:`bool`:
        """

        # Load model
        model = wc_lang.io.Reader().run(model_path)

        # Set concentration of species to 0:
        for tweak_specie_id in tweak_specie_ids:
            temp = re.findall('[a-zA-Z_0-9]+', tweak_specie_id)
            tweak_specie_type_id = temp[0]
            tweak_compartment_id = temp[1]
            tweak_specie_compartment = model.compartments.get_one(id=tweak_compartment_id)
            tweak_specie = model.species_types.get_one(id=tweak_specie_type_id).species.get_one(compartment=tweak_specie_compartment)
            tweak_specie.concentration.value = 0

        # Run model
        results_dir = os.path.expanduser('~/tmp/checkpoints_dir/')
        simulation = Simulation(model)
        num_events, results_dir = simulation.run(end_time = end_time, results_dir = results_dir, checkpoint_period = checkpoint_period)
        run_results = RunResults(results_dir)

        # Check target species concentration at the end of simulation
        is_constant = []
        for target_specie_id in target_specie_ids:
            temp = re.findall('[a-zA-Z_0-9]+', target_specie_id)
            target_specie_type_id = temp[0]
            target_compartment_id = temp[1]
            target_compartment = model.compartments.get_one(id=target_compartment_id)
            target_specie = model.species_types.get_one(id=target_specie_type_id).species.get_one(compartment=target_compartment)
            concentration = run_results.get('populations')[target_specie.id()][:].values #convert panda.series to np.ndarray

            if all(concentration[0]==concentration):
                is_constant.append('True')
            else:
                is_constant.append('False')

        return is_constant, run_results

    def concentrations_scan(self, model_path, end_time, checkpoint_period, tweak_specie_ids, concentrations):
        # Load model
        model = wc_lang.io.Reader().run(model_path)
        run_results =[]

        # Loop through the
        for concentration in concentrations:

            # Set concentration of species to 0:
            for tweak_specie_id in tweak_specie_ids:
                temp = re.findall('[a-zA-Z_0-9]+', tweak_specie_id)
                tweak_specie_type_id = temp[0]
                tweak_compartment_id = temp[1]
                tweak_specie_compartment = model.compartments.get_one(id=tweak_compartment_id)
                tweak_specie = model.species_types.get_one(id=tweak_specie_type_id).species.get_one(compartment=tweak_specie_compartment)
                tweak_specie.concentration.value = concentration

            # Run model
            results_dir = os.path.expanduser('~/tmp/checkpoints_dir/')
            simulation = Simulation(model)
            num_events, results_dir = simulation.run(end_time = end_time, results_dir = results_dir, checkpoint_period = checkpoint_period)
            run_result = RunResults(results_dir)
            run_results.append(run_result)

        return run_results

    def is_constant_tweak_reactions(self, model_path, end_time, checkpoint_period, tweak_reaction_ids, target_specie_ids):
        # Load model
        model = wc_lang.io.Reader().run(model_path)

        # Set concentration of species to 0:
        for tweak_reaction_id in tweak_reaction_ids:
            for reaction in model.get_reactions():
                if reaction.id == tweak_reaction_id:
                    reaction.rate_laws[0].k_m = 0
                    reaction.rate_laws[0].k_cat = 0
                    break

        # Run model
        results_dir = os.path.expanduser('~/tmp/checkpoints_dir/')
        simulation = Simulation(model)
        num_events, results_dir = simulation.run(end_time = end_time,
                                                 results_dir = results_dir,
                                                 checkpoint_period = checkpoint_period)
        run_results = RunResults(results_dir)

        # Check target species concentration at the end of simulation
        is_constant = []
        for target_specie_id in target_specie_ids:
            temp = re.findall('[a-zA-Z_0-9]+', target_specie_id)
            target_specie_type_id = temp[0]
            target_compartment_id = temp[1]
            target_compartment = model.compartments.get_one(id=target_compartment_id)
            target_specie = model.species_types.get_one(id=target_specie_type_id).species.get_one(compartment=target_compartment)
            concentration = run_results.get('populations')[target_specie.id()][:].values #convert panda.series to np.ndarray

            if all(concentration[0]==concentration):
                is_constant.append('True')
            else:
                is_constant.append('False')

        return is_constant, run_results

    def reactions_scan(self, model_path, end_time, checkpoint_period, tweak_reaction_ids, k_cats):
        # Load model
        model = wc_lang.io.Reader().run(model_path)
        run_results=[]

        for k_cat in k_cats:
            # Set concentration of species to 0:
            for tweak_reaction_id in tweak_reaction_ids:
                for reaction in model.get_reactions():
                    if reaction.id == tweak_reaction_id:
                        reaction.rate_laws[0].k_cat = k_cat
                        break

            # Run model
            results_dir = os.path.expanduser('~/tmp/checkpoints_dir/')
            simulation = Simulation(model)
            num_events, results_dir = simulation.run(end_time = end_time,
                                                     results_dir = results_dir,
                                                     checkpoint_period = checkpoint_period)
                                                     
        run_result = RunResults(results_dir)
        run_results.append(run_result)
        return run_results
