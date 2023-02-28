"""Low-level method that interact with the API for SPH computation.

These functions are to be used inside the higher level constructs of the
Inductiva client. For instance, `scenario.simulate()` uses the `run_simulation`
function.
"""
import pathlib
from inductiva.api import invoke_api
from inductiva.types import Path

from typing import Optional

# pylint: disable=unused-argument


def run_simulation(sim_dir: Path,
                   input_file_name: str = "splishsplash_input.json",
                   output_dir: Optional[Path] = None) -> pathlib.Path:
    """Run SplishSplash in the API.

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
