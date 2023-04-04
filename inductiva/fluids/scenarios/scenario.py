"""Base class for scenarios."""

from abc import ABC
from functools import singledispatch, singledispatchmethod
import tempfile
from typing import Optional

from inductiva.types import Path
from inductiva.simulation import Simulator


class Scenario(ABC):
    """Base class for scenarios."""

    def simulate(
        self,
        simulator: Simulator,
        output_dir: Optional[Path] = None,
        **kwargs,
    ):
        """Simulates the scenario."""

        with tempfile.TemporaryDirectory() as input_dir:

            self.gen_aux_files(simulator, input_dir)
            self.gen_config(simulator, input_dir)

            output_path = simulator.run(
                input_dir,
                self.get_config_filename(simulator),
                output_dir=output_dir,
                **kwargs,
            )

        return output_path

    @singledispatchmethod
    @classmethod
    def get_config_filename(cls, simulator: Simulator):
        pass

    @singledispatchmethod
    def gen_aux_files(self, simulator: Simulator, input_dir: str):
        pass

    @singledispatchmethod
    def gen_config(self, simulator: Simulator, input_dir: str):
        pass
