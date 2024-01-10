"""Base class for low-level simulators."""
from typing import Optional
from abc import ABC

from inductiva import types, tasks, resources
from inductiva.utils import files


class Simulator(ABC):
    """Base simulator class."""

    def __init__(self):
        self.api_method_name = ""
        self._is_mpi_available = False

    def override_api_method_prefix(self, prefix: str):
        """Override the API method prefix.

        Example:
            # prefix = "protein_solvation"
            "md.gromacs.run_simulation" becomes
              "protein_solvation.gromacs.run_simulation"

        Args:
            prefix: The new prefix to use.
        """
        last_elements = self.api_method_name.split(".")[1:]
        all_elements = [prefix] + last_elements

        self.api_method_name = ".".join(all_elements)

    def _setup_input_dir(self, input_dir: types.Path):
        """Setup the simulator input directory."""
        input_dir = files.resolve_path(input_dir)
        if not input_dir.is_dir():
            raise ValueError(
                f"The provided path (\"{input_dir}\") is not a directory.")
        return input_dir

    def run(
        self,
        input_dir: types.Path,
        *_args,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[types.Path] = "",
        **kwargs,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory containing the input files.
            _args: Unused in this method, but defined to allow for more
                non-default arguments in method override in subclasses.
            on: The computational resource to launch the simulation in. If None
                the simulation is launched in a machine of the default pool.
            storage_dir: Parent directory for storing simulation
                               results.
            **kwargs: Additional keyword arguments to be passed to the
                simulation API method.
        """
        input_dir = self._setup_input_dir(input_dir)

        if not self._is_mpi_available and isinstance(on, resources.MPICluster):
            raise ValueError("MPI is not available for this simulator. "
                             "Please use a different computational resource.")

        return tasks.run_simulation(
            self.api_method_name,
            input_dir,
            on=on,
            storage_dir=storage_dir,
            **kwargs,
        )
