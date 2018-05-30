""" Methods for verifying models

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-05-10
:Copyright: 2018, Karr Lab
:License: MIT
"""

import os
import wc_lang
import unittest
from wc_sim.multialgorithm.simulation import Simulation
from wc_sim.multialgorithm.run_results import RunResults

class SubmodelDynamicsTestCase(unittest.TestCase):
    """ Methods for verifying submodels  """

    results_dir = os.path.expanduser('~/tmp/checkpoints_dir/')

    def is_species_constant(model_path, end_time, checkpoint_period,
                            tweak_specie_type_id, tweak_specie_compartment_id,
                            target_specie_type_id, target_specie_compartment_id):

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

            TODO: - integrate into test_full_model.py
        """

        # Load model
        model = wc_lang.io.Reader().run(model_path)

        # Tweak concentration fo specie
        tweak_specie_compartment = model.compartments.get_one(id=tweak_specie_compartment_id)
        tweak_specie = model.species_types.get_one(id=tweak_specie_type_id).species.get_one(compartment=tweak_specie_compartment)
        tweak_specie.concentration.value = 0

        # Run model
        simulation = Simulation(model)
        num_events, results_dir = simulation.run(end_time = end_time,
                                                 results_dir = results_dir,
                                                 checkpoint_period = end_time) # Save tiem by having 1 checkpoint
        run_results = RunResults(results_dir)

        # Check target species concentration at the end of simulation
        target_specie_compartment = model.compartments.get_one(id=target_specie_compartment_id)
        target_specie = model.species_types.get_one(id=target_specie_type_id).species.get_one(compartment=target_specie_compartment)
        concentration = run_results.get('populations')[tweak_specie.id()][:].values #convert panda.series to np.ndarray

        if all(concentration[0]==concentration):
            return True
        else:
            return False

    def is_reaction_essential(model, reaction, target_species_type_id, target_compartment_id, time, is_constant):
        pass

    def parameter_scan(model, reaction, target_species_type_id, target_compartment_id, time, is_constant):
        pass
