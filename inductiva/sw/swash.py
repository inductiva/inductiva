"""Low-level method that interacts with the API for SWASH computation."""
import pathlib
from ..api import invoke_api
from ..types import Path

from typing import Optional

# pylint: disable=unused-argument


def run_simulation(sim_dir: Path,
                   input_filename: str,
                   n_cores: int,
                   output_dir: Optional[Path] = None) -> pathlib.Path:
    """Run SWASH simulation in the API.

    Args:
        sim_dir: Path to the directory containing the simulation inputs.

    Returns:
        TODO: once we have a class for holding the outputs of the simulation,
        it should be used here
    """
    params = locals()
    del params["output_dir"]

    return invoke_api(params, run_simulation, output_dir=output_dir)


# pylint: enable=unused-argument
