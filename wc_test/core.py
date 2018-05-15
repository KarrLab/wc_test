""" Methods for verifying models

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-05-10
:Copyright: 2018, Karr Lab
:License: MIT
"""

import unittest


class SubmodelDynamicsTestCase(unittest.TestCase):
    """ Methods for verifying submodels
    """

    def is_species_essential(model, species_type_id, compartment_id, target_species_type_id, target_compartment_id, time, is_constant):
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

        specie = model.species_types.get_one(id=species_type_id).species.get_one(compartment=compartment_id)
        specie.concetration.value = 0

        if is_constant:
            target_concetration = target_specie.concetration.value

        """ add code to run model in wc_sim once Arthur fixed repo """

        target_specie = model.species_types.get_one(id=target_species_type_id).species.get_one(compartment=target_compartment_id)

        if is_constant:
            target_specie.concetration.value
            self.assertAlmostEqual(target_specie.concetration.value, target_concetration)  # results at time t
        else:
            self.assertAlmostEqual(target_specie.concetration.value, 0)

    def is_reaction_essential(model, reaction, target_species_type_id, target_compartment_id, time, is_constant):
        """ Checks whether setting K_cat of 'reaction' to 0 either:

        * reduces the concentration of specie 'target_id' (instance of 'target_class') to 0 within 'time'
        * stabilizes (i.e. remains constant) the concentration of target_species_type[target_compartment] within 'time'

        Args:
            model (:obj:`wc_lang.Model`): model
            reaction (:obj:`wc_lang.Reaction`): reaction
            compartment (:obj:`wc_lang.CompartmentModel`): compartment
            target_species_type (:obj:`wc_lang.SpeciesType`):
            target_compartment (:obj:`wc_lang.CompartmentModel`):
            time (:obj:`float`): time
            is_constant (:obj:`bool`):

        Returns:
            :obj:`bool`:
        """
        pass
