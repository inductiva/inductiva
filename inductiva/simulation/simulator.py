"""Base class for low-level simulators."""
from abc import ABC, abstractmethod
import pathlib
from typing import Any, Dict, Optional, Literal

from inductiva import api
from inductiva import types
from inductiva.utils import files


class Simulator(ABC):
    """Base simulator class."""
    @property
    @abstractmethod
    def api_method_name(self) -> str:
        pass

    def run(self,
            sim_dir: types.Path,
            output_dir: Optional[types.Path] = None,
            **kwargs,
        ) -> pathlib.Path:
        """Run the simulation."""
        sim_dir = pathlib.Path(sim_dir)
        if not sim_dir.is_dir():
            raise ValueError(
                f"The provided path ({sim_dir}) is not a directory.")

        if output_dir is None:
            output_dir = sim_dir.with_name(f"{sim_dir.name}-output")
            output_dir = files.get_timestamped_path(output_dir)
        else:
            output_dir = files.resolve_path(output_dir)

        return api.run_simulation(
            self.api_method_name,
            sim_dir,
            output_dir,
            params=kwargs,
        )

