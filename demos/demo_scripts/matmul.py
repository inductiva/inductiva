"""
Sample usage of the inductiva package.
"""
import inductiva
import numpy as np

from absl import logging
from absl import flags
from absl import app

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")


def main(_):
    inductiva.api_url = FLAGS.api_url

    m = np.random.randint(10, size=(10, 10))
    n = np.random.randint(10, size=(10, 10))

    local_result = np.matmul(m, n)
    remote_result = inductiva.math.matmul(m=m, n=n)

    success = np.allclose(local_result, remote_result)

    logging.info("Operation successful" if success else \
        "Operation unsuccessful")


if __name__ == "__main__":
    app.run(main)
