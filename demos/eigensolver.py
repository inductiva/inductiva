"""Sample usage of SLEPc eigensolver provided by Inductiva's Web API.

This is an example on how to initialize a connection with the Inductiva
Web API and call a function from linalg package to find eigenvalues
and eigenvectors of a given matrix.
"""
from absl import logging

import inductiva
import scipy
import numpy as np
import time

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.update_config(address="http://192.168.1.50:8000", output_dir="output")

    # Creatae a sparse matrix
    start_time = time.perf_counter()
    matrix_size = 1000000
    diags = [
        np.random.normal(size=matrix_size - 1),
        np.random.normal(size=matrix_size),
        np.random.normal(size=matrix_size - 1),
    ]
    matrix = scipy.sparse.diags(diagonals=diags, offsets=[-1, 0, 1])

    remote_result = inductiva.slepc.linalg.eigs(matrix=matrix,
                                                num_eigenpairs=10)
    start_time = time.perf_counter()

    logging.info("Resulting eigenvalues and eigenvectors %s",
                 str(remote_result))
