import copy
import shutil
import tempfile
import unittest
import wc_kb
import wc_lang
from wc_sim.simulation import Simulation
from wc_sim.run_results import RunResults


class KnowledgeBaseTestCase(unittest.TestCase):
    """ Test case for a knowledge base

    Attributes:
        KB (:obj:`wc_kb.KnowledgeBase`): knowledge base
    """
    pass


class ModelTestCase(unittest.TestCase):
    """ Test case for a model

    Class provides

    * Convenience methods for testing models
    * Basic tests applicable to all models

        * Reactions are element and charged balance

    Attributes:
        KB (:obj:`wc_kb.KnowledgeBase`): knowledge base
        MODEL (:obj:`wc_lang.Model`): model
    """

    def test_reactions_balanced(self):
        """ Check that each reaction is element and charge balanced

        Raises:
            :obj:`Exception`: if one or more reactions is unbalanced
        """

        imbalances = []
        for reaction in self.MODEL.reactions:
            imbalance = self.is_reaction_balanced(reaction)
            if imbalance:
                imbalances.append(imbalance)

        if imbalances:
            msg = []
            for imbalance in imbalances:
                if imbalance.element:
                    msg.append('Reaction {} is element imbalanced: {}'.format())
                if imbalance.charge:
                    msg.append('Reaction {} is charge imbalanced: {}'.format())

            raise Exception('The following reactions are not balanced:\n  {}'.format('\n  '.join(msg)))

    def is_reaction_balanced(self, reaction):
        """ Check that a reaction is element and charge balanced

        Returns:
            imbalance
        """
        pass


class SubmodelTestCase(unittest.TestCase):
    """ Test case for a submodel

    Attributes:
        KB (:obj:`wc_kb.KnowledgeBase`): knowledge base
        SUBMODEL (:obj:`wc_lang.Submodel`): submodel
    """
    pass


class ModelSimulationTestCase(unittest.TestCase):
    """ Test case for a simulation of a model

    Provides basic functionality needed to test simulations of models

    * Perturb model parameters
    * Simulate model
    * Log simulation results
    * Read simulation results
    * Reduce simulation results

    Attributes:
        KB (:obj:`wc_kb.KnowledgeBase`): knowledge base
        MODEL (:obj:`wc_lang.Model`): model
        _results_path (:obj:`str`): path to temporarily save simulation result
    """
    def setUp(self):
        self._results_path = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._results_path)

    def simulate(self, end_time, checkpoint_period, peturbation=None):
        if perturbation:
            model = copy.deepcopy(self.MODEL)
            peturbation(model)
        else:
            model = self.MODEL

        sim = Simulation(self.MODEL)
        num_events, _ = sim.run(time_max=end_time,
                                results_dir=self._results_path,
                                checkpoint_period=checkpoint_period)
        return RunResults(self._results_path)

    def simulate_edge_case(self, end_time, checkpoint_period, param_id, value=0):
        def peturbation(model):
            param = model.parameters.get(id=param_id)
            param.value = value
        return self.simulate(end_time, checkpoint_period, peturbation=peturbation)


class SubmodelSimulationTestCase(unittest.TestCase):
    """ Test case for a simulation of a submodel

    Attributes:
        KB (:obj:`wc_kb.KnowledgeBase`): knowledge base
        SUBMODEL (:obj:`wc_lang.Submodel`): submodel
        _results_path (:obj:`str`): path to temporarily save simulation result
    """

    def setUp(self):
        self._results_path = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._results_path)
