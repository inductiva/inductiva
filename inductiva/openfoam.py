"""Low-level method that interacts with the API for OpenFoam computation."""
import pathlib
from .api import invoke_api
from .types import Path

from typing import Optional

# pylint: disable=unused-argument


def run_simulation(sim_dir: Path,
                   openfoam_simulator: str,
                   n_cores: int,
                   output_dir: Optional[Path] = None) -> pathlib.Path:
    """Run OpenFOAM simulation in the API.

    Args:
        sim_dir: Path to the directory containing the simulation inputs.
        openfoam_simulator: 
        n_cores: Number of CPU cores
        openfoam_solver: specific solver to simulate with OpenFOAM. 
            OpenFOAM contains lots of solvers inside of it, which are used
            to call the run simulation through terminal, e.g.,
            [isoFoam, sonicFoam,...] 
    """
    params = locals()
    del params["output_dir"]

    return invoke_api(params, run_simulation, output_dir=output_dir)


# pylint: enable=unused-argument
