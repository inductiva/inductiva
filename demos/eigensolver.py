"""
Sample usage of the inductiva package.
"""
import inductiva
import numpy as np
import scipy

from absl import logging

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.init(address="http://localhost:8000", output_dir="output")

    m = scipy.sparse.random(10, 10, density=0.01)

    remote_result = inductiva.linalg.eigensolver(matrix=m, num_eigenpairs=10)

    logging.info(remote_result)
