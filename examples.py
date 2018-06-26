import pkg_resourcesimport https://urldefense.proofpoint.com/v2/url?u=http-3A__wc-5Fkb.io&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-A&r=k9rpzjRqL_Wxf5JuzSfjOOwxBMxBKkNmbXN4FD-lhss&m=aRtCI0sE7qfYGHvNkV7sYpS7Yxtabfr6Plui-PKc03I&s=0sWxKVpC6uHacj1PVmNT3fcSf6PRGqmEb5eY5xjvJb0&e=import wc_test
# read KBkb_reader = wc_kb.io.Reader()kb = kb_reader.run(    core_path=pkg_resources.resource_filename('mycoplasma_pneumoniae', 'kb', 'core.xlsx'),    seq_path=pkg_resources.resource_filename('mycoplasma_pneumoniae', 'kb', 'seq.fna'))
# generate modelmodel_gen = mycoplasma_pneumoniae.model_gen.core.ModelGenerator(kb)model = model_gen.run()
class KnowledgeBaseTestCase(wc_test.KnowledgeBaseTestCase):    KB = kb
class ModelTestCase(wc_test.ModelTestCase):    KB = kb    MODEL = model
class TranscriptionSubmodelTestCase(wc_test.SubmodelTestCase):    KB = kb    SUBMODEL = model.submodels.get(id='transcription')class ModelSimulationTestCase(wc_test.ModelSimulationTestCase):    KB = kb    MODEL = model
class TranscriptionSubmodelSimulationTestCase(wc_test.SubmodelSimulationTestCase):    KB = kb    SUBMODEL = model.submodels.get(id='translation')