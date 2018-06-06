import wc_test
import numpy as np

model_path = 'tests/fixtures/min_model.xlsx'

init_concentrations = np.linspace(0,20000,6)
print(init_concentrations)
final_concentrations = wc_test.SubmodelDynamicsTestCase().scan_species(model_path, 900, 10, ['RNA_4[c]'], 'RNA_4[c]', init_concentrations)
print(final_concentrations)
