"""Client for Inductiva's web API."""
import os

from . import fluids

from absl import logging

api_url = os.environ.get("INDUCTIVA_API_URL", "http://api.inductiva.ai")
output_dir = os.environ.get("INDUCTIVA_OUTPUT_DIR", "inductiva_output")

logging.set_verbosity(logging.INFO)
