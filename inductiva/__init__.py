"""Client for Inductiva's web API."""
import os
import logging
import urllib3

from . import fluids
from . import core
from . import molecules
from . import tasks
from . import admin

import absl 

api_url = os.environ.get("INDUCTIVA_API_URL", "http://api.inductiva.ai")
output_dir = os.environ.get("INDUCTIVA_OUTPUT_DIR", "inductiva_output")
api_key = os.environ.get("INDUCTIVA_API_KEY")
working_dir = None

absl.logging.set_verbosity(logging.INFO)

urllib3_logger = logging.getLogger("urllib3.connectionpool")
urllib3_logger.setLevel(logging.CRITICAL)
