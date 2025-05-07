"""Base class for low-level simulators."""
from typing import List, Literal, Optional
from abc import ABC
import logging
import os
import re

import pathlib

from inductiva import projects, types, tasks, resources
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

    _supported_device_values = ["gpu", "cpu", "auto"]

    def __init__(self,
                 /,
                 version: Optional[str] = None,
                 use_dev: bool = False,
                 device: Literal["auto", "cpu", "gpu"] = "auto"):
        """Initialize the simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
            device (str): Specifies whether to use the CPU or GPU version of
                the Docker image. If `auto` is picked we will pick the
                appropriate device based on the machine used to run the
                simulation, meaning, if the machine has a GPU we will pick
                GPU, otherwise we pick CPU.
        """
        if version is not None and not isinstance(version, str):
            raise ValueError("Version must be a string or None.")

        if  device is not None and \
            device.lower() not in self._supported_device_values:
            raise ValueError("Wrong value for `device`. Supported "
                             "values are `gpu`, `cpu` or `auto`.")
        self.simulator = ""
        self.simulator_name_alias = None
        self._version = version
        self._use_dev = bool(use_dev)
        self._image_uri = self._get_image_uri()
        self.container_image = self._image_uri
        self._logger.info("")
        self._device = device.lower() if device else None

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

    def _input_files_exist(self, input_dir, remote_assets, **kwargs):
        """
        Checks if all the files in kwargs are present in the input_dir.
        """
        if remote_assets is not None:
            return

        missing_files = []

        for _, file_path in kwargs.items():
            # Get the full file path by joining the input directory with the
            # file path in kwargs
            full_file_path = os.path.join(input_dir, file_path)

            # Check if the file exists
            if not os.path.isfile(full_file_path):
                missing_files.append(file_path)

        if missing_files:
            missing_files_str = ", ".join(missing_files)
            raise FileNotFoundError(
                "The following files are missing from your input directory:\n"
                f"{missing_files_str}")

    def _regex_exists_in_file(self, file_path: str, pattern: str) -> bool:
        """
        Check if a given regular expression exists in a file.
        
        :param file_path: Path to the file.
        :param pattern: Regular expression pattern to search for.
        :return: True if the pattern exists, False otherwise.
        """
        with open(file_path, "r", encoding="utf-8") as file:
            for line in file:
                if re.search(pattern, line):
                    return True
        return False

    def _check_vcpus(
        self,
        n_vcpus: Optional[int],
        machine: types.ComputationalResources,
    ) -> bool:
        """ Checks if the machine supports the number of n_vcpus passed."""
        if n_vcpus is not None and n_vcpus > machine.n_vcpus.total:
            raise ValueError(
                "The number of virtual cpus asked surpasses the"
                " available virtual cpus for the selected resource.")

    def _get_image_uri(self):
        """Get the appropriate image name for this simulator."""

        img_type = "development" if self._use_dev else "production"
        sim_name = self.name
        name = sim_name.lower()

        available = list_available_images()
        listing = available.get(img_type, {}).get(name, [])

        self._supported_versions = {
            element.split("_")[0] for element in listing
        }

        self._supported_versions_with_suffixes = listing

        if self._version is None:
            # kutu does not have a specific tag for the latest version
            # this hack is a workaround to get that version, but it is prone
            # to errors.
            self._version = max(self._supported_versions)

        if self._version not in self._supported_versions:
            raise ValueError(
                f"Version {self.version} is not available for simulator {name}."
                f" Available versions are: {self._supported_versions}.")
        self._logger.info("■ Using %s image of %s version %s", img_type,
                          sim_name, self.version)

        return f"docker://inductiva/kutu:{name}_v{self._version}"

    @classmethod
    def get_supported_resources(cls):
        """Get the supported computational resources for this simulator."""
        return tuple(cls._supported_resources)

    def _setup_input_dir(self, input_dir: str):
        """Setup the simulator input directory."""
        input_dir_path = pathlib.Path(input_dir)
        if not input_dir_path.is_dir():
            raise ValueError(
                f"The provided path (\"{input_dir}\") is not a directory.")
        return input_dir_path

    def _validate_input_files(self, input_dir, remote_assets):
        if input_dir == "":
            raise ValueError(
                "input_dir cannot be an empty string. Use None instead.")

        if not input_dir and not remote_assets:
            raise ValueError(
                "Either input_dir or remote_assets must be provided.")

    def _get_version_suffixes(self,
                              resource: types.ComputationalResources) -> str:
        """
        Gets the suffixes based on the resource used for the simulation.

        If the resource has a GPU we will prioritize that.

        :param resource: The resource used for the simulation.
        :return: A string of suffixes.
        """
        dev_suffix = "_dev" if self._use_dev else ""
        gpu_suffix = "_gpu" if resource.has_gpu() and (
            f"{self.version}_gpu"
            in self._supported_versions_with_suffixes) else ""

        # Overwrites the gpu suffix if the user passed a specific
        # device
        if self._device == "gpu":
            gpu_suffix = "_gpu"
        elif self._device == "cpu":
            gpu_suffix = ""

        # Add new suffixes here
        suffix = f"{gpu_suffix}"

        if self._device and self._device == "gpu" and not resource.has_gpu():
            raise ValueError("The selected machine needs to have a GPU.\n"
                             "If available, either pick device=`cpu` or pick"
                             " a different machine.")

        suffix = f"{suffix}{dev_suffix}"
        return suffix

    def get_simulator_image_based_on_resource(
            self, resource: types.ComputationalResources):
        """
        Get the simulator image for the simulation, based on the resourced used.

        This method takes in consideration the specified device for the
        simulation. If device=`auto` we will pick the appropriate device based
        on the machine.

        :param resource: The computational resource to use for the simulation.
        :return: The simulator image based on the resource used.
        """
        suffixes = ""

        if self._device:
            suffixes = self._get_version_suffixes(resource)

        return f"{self._image_uri}{suffixes}"

    def run(
        self,
        input_dir: Optional[str],
        *_args,
        on: types.ComputationalResources,
        storage_dir: Optional[str] = "",
        resubmit_on_preemption: bool = False,
        remote_assets: Optional[List[str]] = None,
        project: Optional[str] = None,
        **kwargs,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            on: The computational resource to launch the simulation in. If None
                the simulation is launched in a machine of the default pool.
            input_dir: Path to the directory containing the input files.
            _args: Unused in this method, but defined to allow for more
                non-default arguments in method override in subclasses.
            storage_dir: Parent directory for storing simulation
                               results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
            project: Name of the project to which the task will be
                assigned. If None, the task will be assigned to
                the default project. If the project does not exist, it will be
                created.
            **kwargs: Additional keyword arguments to be passed to the
                simulation API method.
        """
        self._validate_input_files(input_dir, remote_assets)

        input_dir_path = self._setup_input_dir(input_dir) if input_dir else None

        if on is None:
            raise ValueError(
                "A default computational resource is no longer "
                "provided. Please specify a computational resource. "
                "Check https://docs.inductiva.ai/en/latest/how_to/"
                "run-parallel_simulations.html "
                "to learn how to create your own computational resource.")

        self.validate_computational_resources(on)

        if "commands" in kwargs:
            cmds = commands.Command.commands_to_dicts(kwargs["commands"])
            kwargs["commands"] = cmds

        # Get the user-specified image name. If not specified,
        # use the default image name for the current simulator
        self._image_uri = kwargs.pop("container_image", self._image_uri)

        # CustomImage does not use suffixes. We can the image as is
        if self.__class__.__name__ != "CustomImage":
            self._image_uri = self.get_simulator_image_based_on_resource(on)

        if on.has_gpu() and "_gpu" not in self._image_uri:
            logging.warning("Attention: The machine you selected has a GPU, but"
                            " the simulator you picked will run on the CPU "
                            "only.\n")

        # This will create the project if it doesn't exist
        if project:
            projects.Project(project)

        return tasks.run_simulation(
            self.simulator,
            input_dir_path,
            simulator_obj=self,
            storage_dir=storage_dir,
            machine_group=on,
            container_image=self._image_uri,
            resubmit_on_preemption=resubmit_on_preemption,
            remote_assets=remote_assets,
            simulator_name_alias=self.simulator_name_alias,
            project_name=project,
            **kwargs,
        )

    def validate_computational_resources(self, resource):
        """Validate the computational resources passed to the run method.

        Args:
            resource: The computational resource to validate.
            valid_resources: The valid computational resources for the simulator
        """

        if not isinstance(resource, tuple(self._supported_resources)):
            raise ValueError(
                "The computational resource is invalid for this simulator. "
                f"Expected one of {self._supported_resources} but got "
                f"{resource}.")
