"""Wind tunnel scenario to run an object over an air flow."""

from functools import singledispatchmethod
import os
import shutil
from typing import Optional, List, Literal
from enum import Enum
from dataclasses import dataclass
from uuid import UUID

from absl import logging

from inductiva.types import Path
from inductiva.scenarios import Scenario
from inductiva.simulation import Simulator
from inductiva.fluids.simulators import OpenFOAM
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       batch_replace_params_in_template)
from inductiva.utils.files import remove_files_with_tag

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "wind_tunnel")
OPENFOAM_TEMPLATE_INPUT_DIR = "openfoam"


@dataclass
class MeshResolution(Enum):
    """Sets particle radius according to resolution."""
    HIGH = [5, 6]
    MEDIUM = [4, 5]
    LOW = [3, 4]


class WindTunnel(Scenario):
    """Physical scenario of a configurable wind tunnel simulation.

    In this scenario, an object is inserted in a wind tunnel described by the
    user. The object is then subject to an air flow determined for which the
    direction and magnitude is defined by the user.
    """

    def __init__(self,
                 flow_velocity: List[float] = None,
                 domain: Optional[dict] = None):
        """Initializes the `WindTunnel` conditions.

        Args:
            flow_velocity (dict): Velocity of the air flow in m/s.
            domain (dict): List containing the lower and upper boundary of
                the wind tunnel in each (x, y, z) direction. It is the
                natural description with the default OpenFOAM simulator.
        """
        if flow_velocity is None:
            logging.info("Using a default flow velocity: [30, 0, 0].")
            self.flow_velocity = [30, 0, 0]
        elif len(flow_velocity) != 3:
            raise ValueError("`flow_velocity` must have 3 values.")
        else:
            self.flow_velocity = flow_velocity

        if domain is None:
            logging.info("Using a default domain: `{\"x\":"
                         "[-5, 15], \"y\": [-4, 4], \"z\": [0, 8]}`.")
            self.domain = {"x": [-5, 15], "y": [-4, 4], "z": [0, 8]}
        elif not isinstance(domain, dict):
            raise ValueError(
                "`domain` must be a dictionary of the type `{\"x\":"
                "[-5, 15], \"y\": [-4, 4], \"z\": [0, 8]}`.")
        else:
            self.domain = domain

    def simulate(self,
                 simulator: Simulator = OpenFOAM(
                     "windtunnel.openfoam.run_simulation"),
                 output_dir: Optional[Path] = None,
                 resource_pool_id: Optional[UUID] = None,
                 object_path: Optional[Path] = None,
                 simulation_time: float = 100,
                 output_time_step: float = 50,
                 resolution: Literal["high", "medium", "low"] = "medium",
                 n_cores: int = 1):
        """Simulates the wind tunnel scenario synchronously.

        Args:
            object_path: Path to object inserted in the wind tunnel.
            simulator: Simulator to use for the simulation.
            output_dir: Path to the directory where the simulation output
                is downloaded.
            simulation_time: Simulation time, in seconds.
            write_interval: Interval between simulation outputs, in seconds.
            n_cores: Number of cores to use for the simulation.
            """

        if object_path:
            self.object_path = object_path
        else:
            logging.info("WindTunnel is empty. Object path not specified.")

        self.simulation_time = simulation_time
        self.output_time_step = output_time_step
        self.n_cores = n_cores
        self.resolution = MeshResolution[resolution.upper()].value

        output_path = super().simulate(
            simulator,
            output_dir=output_dir,
            resource_pool_id=resource_pool_id,
            n_cores=n_cores,
            method_name="simpleFoam",
        )

        return output_path

    def simulate_async(self,
                       simulator: Simulator = OpenFOAM(
                           "windtunnel.openfoam.run_simulation"),
                       resource_pool_id: Optional[UUID] = None,
                       object_path: Optional[Path] = None,
                       simulation_time: float = 100,
                       output_time_step: float = 50,
                       resolution: Literal["high", "medium", "low"] = "medium",
                       n_cores: int = 1):
        """Simulates the wind tunnel scenario asynchronously.

        Args:
            object_path: Path to object inserted in the wind tunnel.
            simulator: Simulator to use for the simulation.
            output_dir: Path to the directory where the simulation output
                is downloaded.
            simulation_time: Simulation time, in seconds.
            write_interval: Interval between simulation outputs, in seconds.
            n_cores: Number of cores to use for the simulation.
            """

        if object_path:
            self.object_path = object_path
        else:
            logging.info("WindTunnel is empty. Object path not specified.")

        self.resolution = MeshResolution[resolution.upper()].value
        self.simulation_time = simulation_time
        self.output_time_step = output_time_step
        self.n_cores = n_cores

        task_id = super().simulate_async(
            simulator,
            resource_pool_id=resource_pool_id,
            n_cores=n_cores,
            method_name="simpleFoam",
        )

        return task_id

    @singledispatchmethod
    @classmethod
    def get_config_filename(cls, simulator: Simulator):  # pylint: disable=unused-argument
        raise ValueError(
            f"Simulator not supported for `{cls.__name__}` scenario.")

    @singledispatchmethod
    def gen_aux_files(self, simulator: Simulator, input_dir: str):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )

    @singledispatchmethod
    def gen_config(self, simulator: Simulator, input_dir: str):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )


@WindTunnel.get_config_filename.register
def _(self, simulator: OpenFOAM):  # pylint: disable=unused-argument
    pass


@WindTunnel.gen_aux_files.register
def _(self, simulator: OpenFOAM, input_dir):  # pylint: disable=unused-argument
    # The WindTunnel with OpenFOAM requires changing multiple files
    template_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                OPENFOAM_TEMPLATE_INPUT_DIR)

    # Copy all files from the template dir to the input directory
    for directory in os.listdir(template_dir):
        shutil.copytree(os.path.join(template_dir, directory),
                        os.path.join(input_dir, directory))

    # Remove all files that have .jinja in input_dir
    remove_files_with_tag(input_dir, ".jinja")


@WindTunnel.gen_config.register
def _(self, simulator: OpenFOAM, input_dir: str):  # pylint: disable=unused-argument
    """Generates the configuration files for OpenFOAM."""

    # The WindTunnel with OpenFOAM requires changing multiple files
    template_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                OPENFOAM_TEMPLATE_INPUT_DIR)

    # Add object path to its respective place
    object_input_dir_path = os.path.join(input_dir, "constant", "triSurface")
    os.mkdir(object_input_dir_path)
    shutil.copy(self.object_path,
                os.path.join(object_input_dir_path, "object.obj"))

    batch_replace_params_in_template(
        templates_dir=template_dir,
        template_filename_paths=[
            os.path.join("0", "include",
                         "initialConditions_template.openfoam.jinja"),
            os.path.join("system", "controlDict_template.openfoam.jinja"),
            os.path.join("system", "blockMeshDict_template.openfoam.jinja"),
            os.path.join("system", "decomposeParDict_template.openfoam.jinja"),
            os.path.join("system", "snappyHexMeshDict_template.openfoam.jinja")
        ],
        params={
            "flow_velocity": self.flow_velocity,
            "simulation_time": self.simulation_time,
            "output_time_step": self.output_time_step,
            "n_cores": self.n_cores,
            "domain": self.domain,
            "resolution": self.resolution,
        },
        output_filename_paths=[
            os.path.join(input_dir, "0", "include", "initialConditions"),
            os.path.join(input_dir, "system", "controlDict"),
            os.path.join(input_dir, "system", "blockMeshDict"),
            os.path.join(input_dir, "system", "decomposeParDict"),
            os.path.join(input_dir, "system", "snappyHexMeshDict")
        ])
