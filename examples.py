import pkg_resourcesimport wc_kbimport wc_kb.ioimport wc_test# read KBkb_reader = wc_kb.io.Reader()kb = kb_reader.run(    pkg_resources.resource_filename('mycoplasma_pneumoniae', 'kb', 'core.xlsx'),    seq_path=pkg_resources.resource_filename('mycoplasma_pneumoniae', 'kb', 'seq.fna'))[wc_kb.KnowledgeBase][0]# generate modelmodel_gen = mycoplasma_pneumoniae.model_gen.core.ModelGenerator(kb)model = model_gen.run()class KnowledgeBaseTestCase(wc_test.KnowledgeBaseTestCase):    KB = kbclass ModelTestCase(wc_test.ModelTestCase):    KB = kb    MODEL = modelclass TranscriptionSubmodelTestCase(wc_test.SubmodelTestCase):    KB = kb    SUBMODEL = model.submodels.get(id='transcription')class ModelSimulationTestCase(wc_test.ModelSimulationTestCase):    KB = kb    MODEL = modelclass TranscriptionSubmodelSimulationTestCase(wc_test.SubmodelSimulationTestCase):    KB = kb    SUBMODEL = model.submodels.get(id='translation')