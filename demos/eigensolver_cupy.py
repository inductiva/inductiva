"""Sample usage of the eigensolver method using CuPy.

This is an example on how to initialize a connection the Inductiva
Web API and call a function from linalg package to find eigenvalues
and eigenvectors of a given matrix.
"""
import inductiva
import scipy
import numpy as np
import time
import utils

from absl import logging
from absl import app
from absl import flags

FLAGS = flags.FLAGS

flags.DEFINE_integer("size", 1000, "Size of the square matrix to use.")

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")


def main(_):
    inductiva.api_url = FLAGS.api_url

    m = utils.get_square_tridiagonal_h_matrix(FLAGS.size)

    time_start = time.perf_counter()
    remote_result = inductiva.cupy.linalg.eigs(m=m)
    logging.info("API time: %s", time.perf_counter() - time_start)

    time_start = time.perf_counter()
    local_result = scipy.sparse.linalg.eigsh(m)
    logging.info("Local time (scipy): %s", time.perf_counter() - time_start)

    logging.info(np.allclose(remote_result[0], local_result[0]))


if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)
    app.run(main)
