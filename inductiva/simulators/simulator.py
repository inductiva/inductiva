"""Base class for low-level simulators."""
from typing import Optional
from abc import ABC
import logging

import pathlib

from inductiva import types, tasks, resources
from .methods import list_available_images
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
    _logger = logging.getLogger(__name__)

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        if version is not None and not isinstance(version, str):
            raise ValueError("Version must be a string or None.")
        self.api_method_name = ""
        self._version = version
        self._use_dev = bool(use_dev)
        self._image_uri = self._get_image_uri()
        self._logger.info("")

    @property
    def version(self):
        """Get the version of the simulator."""
        return self._version

    @property
    def use_dev(self):
        """Get whether the development version of the simulator is used."""
        return self._use_dev

    @property
    def name(self):
        """Get the name of the simulator."""
        return self.__class__.__name__

    @property
    def image_uri(self):
        """Get the image URI for this simulator."""
        return self._image_uri

    def _get_image_uri(self):
        """Get the appropriate image name for this simulator."""

        img_type = "development" if self._use_dev else "production"
        sim_name = self.name
        name = sim_name.lower()

        available = list_available_images()
        listing = available.get(img_type, {}).get(name, [])

        if self._version is None:
            # kutu does not have a specific tag for the latest version
            # this hack is a workaround to get that version, but it is prone
            # to errors.
            self._version = max(listing)

        if self._version not in listing:
            raise ValueError(
                f"Version {self.version} is not available for simulator {name}."
                f" Available versions are: {listing}.")
        self._logger.info("■ Using %s image of %s version %s", img_type,
                          sim_name, self.version)

        suffix = "_dev" if self._use_dev else ""
        return f"docker://inductiva/kutu:{name}_v{self._version}" + suffix

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

    def _setup_input_dir(self, input_dir: str):
        """Setup the simulator input directory."""
        input_dir_path = pathlib.Path(input_dir)
        if not input_dir_path.is_dir():
            raise ValueError(
                f"The provided path (\"{input_dir}\") is not a directory.")
        return input_dir_path

    def run(
        self,
        input_dir: str,
        *_args,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[str] = "",
        resubmit_on_preemption: bool = False,
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
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiates with
                `spot=True`.
            **kwargs: Additional keyword arguments to be passed to the
                simulation API method.
        """
        input_dir_path = self._setup_input_dir(input_dir)

        self.validate_computational_resources(on)

        if "commands" in kwargs:
            cmds = commands.Command.commands_to_dicts(kwargs["commands"])
            kwargs["commands"] = cmds

        # Get the user-specified image name. If not specified,
        # use the default image name for the current simulator
        container_image = kwargs.pop("container_image", self._image_uri)

        return tasks.run_simulation(
            self.api_method_name,
            input_dir_path,
            simulator=self,
            storage_dir=storage_dir,
            computational_resources=on,
            extra_metadata=extra_metadata,
            container_image=container_image,
            resubmit_on_preemption=resubmit_on_preemption,
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
