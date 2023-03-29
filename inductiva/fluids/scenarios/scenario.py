"""Base class for scenarios."""

from abc import ABC
from typing import Optional

from inductiva.types import Path
from inductiva.simulation import Simulator


class Scenario(ABC):
    """Base class for scenarios."""

    def simulate(
        self,
        simulator: Simulator,
        sim_dir: Path,
        *args,
        output_dir: Optional[Path] = None,
        **kwargs,
    ):
        """Simulates the scenario."""

        output_path = simulator.run(sim_dir,
                                    args,
                                    output_dir=output_dir,
                                    kwargs=kwargs)

        return output_path
