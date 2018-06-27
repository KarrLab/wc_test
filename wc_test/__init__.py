import pkg_resources
from .core import ModelTestCase
from .core import KnowledgeBaseTestCase
from .core import StaticTestCase
from .core import DynamicTestCase

# read version
with open(pkg_resources.resource_filename('wc_test', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
