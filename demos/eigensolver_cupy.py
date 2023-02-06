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

    inductiva.init(address="http://localhost:8080", output_dir="output")

    logging.info("size ...")

    size = 10000000
    logging.info("random ...")

    diags = [
        np.random.normal(size=size - 1),
        np.random.normal(size=size),
        np.random.normal(size=size - 1),
    ]

    m = scipy.sparse.diags(diagonals=diags, offsets=[-1, 0, 1], format="csr")

    logging.info("sum ...")

    m = m + m.transpose()

    logging.info("Starting ...")

    time_start = time.perf_counter()
    remote_result = inductiva.linalg.eigs_cupy(m=m)
    logging.info(time.perf_counter() - time_start)


    time_start = time.perf_counter()
    local_result = scipy.sparse.linalg.eigsh(m)
    logging.info(time.perf_counter() - time_start)

    logging.info(np.allclose(remote_result[0], local_result[0]))
