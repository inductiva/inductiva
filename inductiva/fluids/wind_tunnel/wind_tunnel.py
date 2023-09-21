"""Wind tunnel scenario to run an object over an air flow."""

from dataclasses import dataclass
import enum
from functools import singledispatchmethod
import os
import shutil
import io
from typing import Optional, List, Literal

from absl import logging
import numpy as np

from inductiva import tasks, resources, types, simulators, scenarios, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "wind_tunnel")
OPENFOAM_TEMPLATE_INPUT_DIR = "openfoam"
FILES_SUBDIR = "files"
COMMANDS_TEMPLATE_FILE_NAME = "commands.json.jinja"


@dataclass
class MeshResolution(enum.Enum):
    """Sets particle radius according to resolution."""
    HIGH = [5, 6]
    MEDIUM = [4, 5]
    LOW = [3, 4]
    VERY_LOW = [2, 3]


class WindTunnel(scenarios.Scenario):
    """Physical scenario of a configurable wind tunnel simulation.

    A wind tunnel is a tool used in aerodynamic research to study the
    effects of air moving past solid objects. Here, the tunnel consists
    of a box object in 3D space (x, y, z) space, where air flows in the
    positive x-direction with a certain velocity.

    An arbitrary object is placed within the tunnel, sucht that air flows
    around it, as illustrated in the schematic below:
    |--------------------------------|
    |->          _____               |
    |->        _/     |              |
    |->_______|_o___O_|______________|

    This scenario solves steady-state continuity and momentum equations
    (time-independent) with incompressible flow.
    The simulation solves the time-independent equations for several
    time steps, based on the state of the previous one. The end goal is
    to determine the steady-state of the system, i.e., where the flow
    does not change in time anymore.

    Currently, the following variables are fixed:
    - The fluid being inject is air.
    - The flow is incompressible (this restricts the max air velocity).
    - Air only flows in the positive x-direction.
    - Some post-processing of the data occurs at run-time: streamlines,
    pressure_field, cutting planes and force coefficients.
    """

    valid_simulators = [simulators.OpenFOAM]

    def __init__(self,
                 flow_velocity: List[float] = None,
                 domain: Optional[dict] = None):
        """Initializes the `WindTunnel` conditions.

        Args:
            flow_velocity (List): Velocity of the air flow (m/s).
                The maximum velocity magnitude allowed is 100 m/s
            domain (dict): List containing the lower and upper boundary of
                the wind tunnel in each (x, y, z) direction (m). It is the
                natural description with the default OpenFOAM simulator.
        """
        if flow_velocity is None:
            logging.info("Using a default flow velocity: [30, 0, 0].")
            self.flow_velocity = [30, 0, 0]
        elif len(flow_velocity) != 3:
            raise ValueError("`flow_velocity` must have 3 values.")
        elif np.linalg.norm(np.array(flow_velocity)) > 100:
            raise ValueError("The `flow_velocity` magnitude is too high,"
                             " it must be less than 100 m/s.")
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
                 simulator: simulators.Simulator = simulators.OpenFOAM(),
                 machine_group: Optional[resources.MachineGroup] = None,
                 run_async: bool = False,
                 object_path: Optional[types.Path] = None,
                 num_iterations: float = 100,
                 resolution: Literal["high", "medium", "low"] = "medium",
                 n_cores: int = 2) -> tasks.Task:
        """Simulates the wind tunnel scenario synchronously.

        Args:
            simulator: Simulator used to simulate the scenario.
                Valid simulators: OpenFOAM.
            object_path: Path to object inserted in the wind tunnel.
            output_dir: Path to the directory where the simulation output
                is downloaded when running synchronously.
            num_iterations: Number of iterations of the simulator.
            n_cores: Number of cores to use for the simulation.
            resolution: Level of detail of the mesh used for the simulation.
                Options: "high", "medium" or "low".
            machine_group: The machine group to use for the simulation.
        """
        simulator.override_api_method_prefix("wind_tunnel")

        if object_path is None:
            raise ValueError("WindTunnel is empty. Object path not specified.")

        self.object_path = utils.files.resolve_path(object_path)
        self.num_iterations = num_iterations
        self.n_cores = n_cores
        self.resolution = MeshResolution[resolution.upper()].value

        commands = self.get_commands()

        task = super().simulate(simulator,
                                machine_group=machine_group,
                                run_async=run_async,
                                n_cores=n_cores,
                                commands=commands)

        return task

    def get_commands(self):
        """Returns the commands for the simulation."""

        commands_template_path = os.path.join(SCENARIO_TEMPLATE_DIR,
                                              OPENFOAM_TEMPLATE_INPUT_DIR,
                                              COMMANDS_TEMPLATE_FILE_NAME)

        inmemory_file = io.StringIO()
        utils.templates.replace_params(
            template_path=commands_template_path,
            params={"n_cores": self.n_cores},
            output_file=inmemory_file,
        )
        commands = self.read_commands_from_file(inmemory_file)

        return commands

    @singledispatchmethod
    def create_input_files(self, simulator: simulators.Simulator):
        pass


@WindTunnel.create_input_files.register
def _(self, simulator: simulators.OpenFOAM, input_dir):  # pylint: disable=unused-argument
    """Creates OpenFOAM simulation input files."""

    # The WindTunnel with OpenFOAM requires changing multiple files
    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      OPENFOAM_TEMPLATE_INPUT_DIR, FILES_SUBDIR)

    # Copy all files from the template dir to the input directory
    shutil.copytree(template_files_dir,
                    input_dir,
                    dirs_exist_ok=True,
                    symlinks=True)

    utils.templates.batch_replace_params(
        templates_dir=input_dir,
        template_filenames=[
            os.path.join("0", "include",
                         "initialConditions_template.openfoam.jinja"),
            os.path.join("system", "controlDict_template.openfoam.jinja"),
            os.path.join("system", "blockMeshDict_template.openfoam.jinja"),
            os.path.join("system", "decomposeParDict_template.openfoam.jinja"),
            os.path.join("system", "snappyHexMeshDict_template.openfoam.jinja")
        ],
        params={
            "flow_velocity": self.flow_velocity,
            "num_iterations": self.num_iterations,
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
        ],
        remove_templates=True,
    )

    # Add object path to its respective place
    object_dir = os.path.join(input_dir, "constant", "triSurface")
    os.mkdir(object_dir)
    shutil.copy(self.object_path, os.path.join(object_dir, "object.obj"))
