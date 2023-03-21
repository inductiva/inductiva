"""Low-level method that interacts with the API for XBeach computation."""
import pathlib
from ..api import invoke_api
from ..types import Path

from typing import Optional


# pylint: disable=unused-argument
def run_simulation(sim_dir: Path,
                   input_filename: str,
                   n_cores: int,
                   output_dir: Optional[Path] = None) -> pathlib.Path:
    """Run XBeach simulation in the API.

    Args:
        sim_dir: Path to the directory containing the simulation inputs.

    Returns:
        TODO(https://github.com/inductiva/inductiva/issues/99): Add class
            for holding the outputs of the simulation.
    """
    params = locals()
    del params["output_dir"]

    return invoke_api(params, run_simulation, output_dir=output_dir)
