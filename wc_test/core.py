""" wc_test

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-5-10
:Copyright: 2018, Karr Lab
:License: MIT
"""

class SubmodelDynamicsTesting(unittest.TestCase):
    """ Class provides environment and methods to test the dynamics of submodels

    Attributes:
        attr_1 (:obj:`type of attr_1`): description of attr_1
        attr_2 (:obj:`type of attr_2`): description of attr_2
        ...
    """

    def __init__(self, arg_1, arg_2, kwarg_1=None, kwarg_2=None):
        """
        Args:
            arg_1 (:obj:`type of arg_1`): description of arg_1
            arg_2 (:obj:`type of arg_2`): description of arg_2
            kwarg_1 (:obj:`type of kwarg_1`, optional): description of kwarg_1
            kwarg_2 (:obj:`type of kwarg_2`, optional): description of kwarg_2
            ...
        """
        self.attr_1 = arg_1
        self.attr_2 = arg_2

    def IsSpecieEssential(model, species_type_id, compartment_id, target_species_type_id, target_compartment_id, time, isconstant):
        """ Checks whether setting the concentration of species_type[compartment] to 0 either:
            - reduces the concentration of specie 'target_id' (instance of 'target_class') to 0 within 'time'
            - stabilizes (i.e. remains constant) the concentration of target_species_type[target_compartment] within 'time'

        Args:
            model (:obj:`wc_lang.Model`):
            species_type (:obj:`wc_lang.SpeciesType`):
            compartment (:obj:`wc_lang.CompartmentModel`):
            target_species_type (:obj:`wc_lang.SpeciesType`):
            target_compartment (:obj:`wc_lang.CompartmentModel`):
            time (:obj:`float`):
            isconstant (:obj:`boolean`):
        Returns:
            :obj:`NoneType`:
        Raises:
            :obj:`NoneType`:
        """

        specie = model.species_types.get_one(id=species_type_id).species.get_one(compartment=compartment_id)
        specie.concetration.value = 0

        if isconstant:
            target_concetration = target_specie.concetration.value

        """ add code to run model in wc_sim once Arthur fixed repo """

        target_specie = model.species_types.get_one(id=target_species_type_id).species.get_one(compartment=target_compartment_id)

        if isconstant:
            target_specie.concetration.value
            self.assertAlmostEqual(target_specie.concetration.value, target_concetration) #results at time t
        else:
            self.assertAlmostEqual(target_specie.concetration.value, 0)

    def IsReactionEssential(model, reaction, target_species_type_id, target_compartment_id, time, isconstant):
        """ Checks whether setting K_cat of 'reaction' to 0 either:
            - reduces the concentration of specie 'target_id' (instance of 'target_class') to 0 within 'time'
            - stabilizes (i.e. remains constant) the concentration of target_species_type[target_compartment] within 'time'

        Args:
            model (:obj:`wc_lang.Model`):
            reaction (:obj:`wc_lang.Reaction`):
            compartment (:obj:`wc_lang.CompartmentModel`):
            target_species_type (:obj:`wc_lang.SpeciesType`):
            target_compartment (:obj:`wc_lang.CompartmentModel`):
            time (:obj:`float`):
            isconstant (:obj:`boolean`):
        Returns:
            :obj:`NoneType`:
        Raises:
            :obj:`NoneType`:
        """
