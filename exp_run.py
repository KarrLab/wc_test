import wc_test
import wc_lang

model_path = 'tests/fixtures/min_model.xlsx'

mass_balanced   = wc_test.SubmodelDynamicsTestCase().are_reactions_mass_balanced(model_path)
charge_balanced = wc_test.SubmodelDynamicsTestCase().are_reactions_charge_balanced(model_path)
