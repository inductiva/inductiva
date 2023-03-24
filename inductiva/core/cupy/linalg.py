"""Linear algebra methods that are computed using CuPy via the API."""
import numpy as np
import scipy
import scipy.sparse
from ..api import invoke_api

from typing import Tuple
# pylint: disable=unused-argument


def eigs(m: scipy.sparse) -> Tuple[np.ndarray, np.ndarray]:
    return invoke_api(locals(), eigs)


# pylint: enable=unused-argument
