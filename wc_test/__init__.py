import pkg_resources
from .core import KnowledgeBaseTestCase
from .core import ModelGenerationTestCase
from .core import ModelStaticTestCase
from .core import ModelDynamicsTestCase

# read version
with open(pkg_resources.resource_filename('wc_test', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
