"""Sample usage of the inductiva linalg package.

This is an example on how to initialize a connection the Inductiva
Web API and call a function from linalg package to find eigenvalues
and eigenvectors of a given matrix.
"""
import inductiva
import scipy

from absl import logging

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.init(address="http://localhost:8000", output_dir="output")

    m = scipy.sparse.random(10, 10, density=0.01)

    remote_result = inductiva.linalg.eigs(matrix=m, num_eigenpairs=10)

    logging.info("Resulting eigenvalues and eigenvectors %s", str(remote_result))
