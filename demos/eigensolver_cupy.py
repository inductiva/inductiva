"""Sample usage of the inductiva linalg package.

This is an example on how to initialize a connection the Inductiva
Web API and call a function from linalg package to find eigenvalues
and eigenvectors of a given matrix.
"""
import inductiva
import scipy
import numpy as np
import time

from absl import logging

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.init(address="http://localhost:8000", output_dir="output")

    size = 1000

    diags = [
        np.random.normal(size=size - 1),
        np.random.normal(size=size),
        np.random.normal(size=size - 1),
    ]

    m = scipy.sparse.diags(diagonals=diags, offsets=[-1, 0, 1], format="csr")

    m = m + m.transpose()

    time_start = time.perf_counter()
    remote_result = inductiva.cupy.linalg.eigs(m=m)
    logging.info("API time: %s", time.perf_counter() - time_start)

    time_start = time.perf_counter()
    local_result = scipy.sparse.linalg.eigsh(m)
    logging.info("Local time (scipy): %s", time.perf_counter() - time_start)

    logging.info(np.allclose(remote_result[0], local_result[0]))
