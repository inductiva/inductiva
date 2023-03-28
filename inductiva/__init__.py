"""Client for Inductiva's web API."""
import os

from . import fluids
from . import core

from absl import logging

api_url = os.environ.get("INDUCTIVA_API_URL", "http://api.inductiva.ai")
output_dir = os.environ.get("INDUCTIVA_OUTPUT_DIR", "inductiva_output")
api_key = os.environ.get("INDUCTIVA_API_KEY")

logging.set_verbosity(logging.INFO)
