"""Base class for low-level simulators."""
from typing import Optional
from abc import ABC

import pathlib

from inductiva import types, tasks, resources
from inductiva import commands


def mpi_enabled(cls):
    """Class decorator that adds MPICluster to the supported resources.
    """
    supports = getattr(cls, "_supported_resources", set())
    cls._supported_resources = supports | {resources.MPICluster}  # pylint: disable=protected-access

    return cls


class Simulator(ABC):
    """Base simulator class."""

    _supported_resources = {
        resources.MachineGroup, resources.ElasticMachineGroup
    }

    def __init__(self):
        self.api_method_name = ""

    @classmethod
    def get_supported_resources(cls):
        """Get the supported computational resources for this simulator."""
        return tuple(cls._supported_resources)

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
        input_dir = pathlib.Path(input_dir)
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
        extra_metadata: Optional[dict] = None,
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

        self.validate_computational_resources(on)

        if "commands" in kwargs:
            kwargs["commands"] = commands.Command.commands_to_dicts(
                kwargs["commands"])

        return tasks.run_simulation(
            self.api_method_name,
            input_dir,
            computational_resources=on,
            storage_dir=storage_dir,
            extra_metadata=extra_metadata,
            **kwargs,
        )

    def validate_computational_resources(self, resource):
        """Validate the computational resources passed to the run method.

        Args:
            resource: The computational resource to validate.
            valid_resources: The valid computational resources for the simulator
        """

        if resource is not None \
            and not isinstance(resource, tuple(self._supported_resources)):
            raise ValueError(
                "The computational resource is invalid for this simulator. "
                f"Expected one of {self._supported_resources} but got "
                f"{resource}.")
