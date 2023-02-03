"""
Sample usage of the inductiva package.
"""
import inductiva
import numpy as np

from absl import logging

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.init(address="http://192.168.1.50:8000", output_dir="output")

    m = np.random.randint(10, size=(10, 10))
    n = np.random.randint(10, size=(10, 10))

    local_result = np.matmul(m, n)
    remote_result = inductiva.math.matmul(m=m, n=n)

    success = np.allclose(local_result, remote_result)

    logging.debug("Operation successful" if success else \
        "Operation unsuccessful")
