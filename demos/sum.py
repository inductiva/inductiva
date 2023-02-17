"""
Sample usage of the inductiva package.
"""
import inductiva

from absl import logging
from absl import flags

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.update_config(FLAGS.api_url)

    result = inductiva.math.sum(a=1, b=1)

    logging.debug("Result = %s", result)
