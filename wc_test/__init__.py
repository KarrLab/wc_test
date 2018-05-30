import pkg_resources
from .core import SubmodelDynamicsTestCase

# read version
with open(pkg_resources.resource_filename('wc_test', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
