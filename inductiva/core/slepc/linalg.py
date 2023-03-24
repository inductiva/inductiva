"""Functions that use the API to compute eigenvalues.

These functions compute the eigenvalues and eigenvectors of
sparse metrices with the corresponding solvers.
"""

from typing import Tuple

import numpy as np
import scipy
from api import invoke_api_f
# pylint: disable=unused-argument


def eigs(matrix: scipy.sparse,
         num_eigenpairs: int) -> Tuple[np.ndarray, np.ndarray]:
    return invoke_api(locals(), eigs)


# pylint: enable=unused-argument
