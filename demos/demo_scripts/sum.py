"""
Sample usage of the inductiva package.
"""
import inductiva

from absl import logging
from absl import flags
from absl import app

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")


def main(_):
    inductiva.api_url = FLAGS.api_url

    result = inductiva.math.sum(a=1, b=1)

    logging.info("Result = %s", result)


if __name__ == "__main__":
    app.run(main)
