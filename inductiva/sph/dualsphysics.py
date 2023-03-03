"""Low-level method that interacts with the API for SPH computation.

These functions are to be used inside the higher level constructs of the
Inductiva client. For instance, `scenario.simulate()` uses the `run_simulation`
function.
"""
import pathlib
from typing import Optional

from inductiva.api import invoke_api
from inductiva.types import Path

# pylint: disable=unused-argument


def run_simulation(sim_dir: Path,
                   input_filename: str,
                   device: str,
                   output_dir: Optional[Path] = None) -> pathlib.Path:
    """Run DualSPHysics in the API.

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
