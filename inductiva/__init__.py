"""Client for Inductiva's web API."""
import os
import logging

import absl

from . import admin
from . import api
from . import core
from . import fluids
from . import generative
from . import molecules
from . import resources
from . import simulation
from . import structures
from . import tasks
from . import templates
from . import utils
from . import stellarators
from . import world
from inductiva.tasks import Task

api_url = os.environ.get("INDUCTIVA_API_URL", "http://api.inductiva.ai")
output_dir = os.environ.get("INDUCTIVA_OUTPUT_DIR", "inductiva_output")
api_key = os.environ.get("INDUCTIVA_API_KEY")
working_dir = None

absl.logging.set_verbosity(logging.INFO)

urllib3_logger = logging.getLogger("urllib3.connectionpool")
urllib3_logger.setLevel(logging.CRITICAL)
