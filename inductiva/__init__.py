"""Client for Inductiva's web API."""
import os
import logging

import absl

from . import api
from . import output_consts
from . import simulators
from . import fluids
from . import coastal
from . import scenarios
from . import generative
from . import molecules
from . import resources
# from . import structures
from . import templates
from . import utils
from . import plasmas
from . import world
from . import tasks

api_url = os.environ.get("INDUCTIVA_API_URL", "http://api.inductiva.ai")
output_dir = os.environ.get("INDUCTIVA_OUTPUT_DIR", "inductiva_output")
api_key = os.environ.get("INDUCTIVA_API_KEY")
working_dir = None

# Disable urllib3 warnings.
# TODO: Verify and fix the appearance of this warning.
urllib3_logger = logging.getLogger("urllib3.connectionpool")
urllib3_logger.setLevel(logging.CRITICAL)

absl.logging.set_verbosity(absl.logging.INFO)

__version__ = "0.3.6"
