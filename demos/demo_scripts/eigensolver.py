"""Sample usage of SLEPc eigensolver provided by Inductiva's Web API.

This is an example on how to initialize a connection with the Inductiva
Web API and call a function from linalg package to find eigenvalues
and eigenvectors of a given matrix.
"""
import scipy
import numpy as np
import time

from absl import logging
from absl import flags
from absl import app

import inductiva
from . import utils

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")


def main(_):
    inductiva.api_url = FLAGS.api_url

    m = utils.get_square_tridiagonal_h_matrix(10)

    time_start = time.perf_counter()
    remote_result = inductiva.slepc.linalg.eigs(matrix=m, num_eigenpairs=10)
    logging.info("API time: %s", time.perf_counter() - time_start)

    time_start = time.perf_counter()
    local_result = scipy.sparse.linalg.eigsh(m)
    logging.info("Local time (scipy): %s", time.perf_counter() - time_start)

    logging.info(np.allclose(remote_result[0], local_result[0]))


if __name__ == "__main__":
    app.run(main)
