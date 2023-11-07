"""Wind tunnel scenario to run an object over an air flow."""

from dataclasses import dataclass
import enum
import os
import shutil
from typing import Optional, List, Literal

from absl import logging
import numpy as np

from inductiva import tasks, resources, types, simulators, scenarios, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "wind_tunnel")
OPENFOAM_TEMPLATE_INPUT_DIR = "openfoam"


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
    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      OPENFOAM_TEMPLATE_INPUT_DIR)

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
            self.params["flow_velocity"] = [30, 0, 0]
        elif len(flow_velocity) != 3:
            raise ValueError("`flow_velocity` must have 3 values.")
        elif np.linalg.norm(np.array(flow_velocity)) > 100:
            raise ValueError("The `flow_velocity` magnitude is too high,"
                             " it must be less than 100 m/s.")
        else:
            self.params["flow_velocity"] = flow_velocity

        if domain is None:
            logging.info("Using a default domain: `{\"x\":"
                         "[-5, 15], \"y\": [-4, 4], \"z\": [0, 8]}`.")
            self.params["domain"] = {"x": [-5, 15], "y": [-4, 4], "z": [0, 8]}
        elif not isinstance(domain, dict):
            raise ValueError(
                "`domain` must be a dictionary of the type `{\"x\":"
                "[-5, 15], \"y\": [-4, 4], \"z\": [0, 8]}`.")
        else:
            self.params["domain"] = domain

    def simulate(
            self,
            simulator: simulators.Simulator = simulators.OpenFOAM(),
            machine_group: Optional[resources.MachineGroup] = None,
            storage_dir: Optional[str] = "",
            object_path: Optional[types.Path] = None,
            num_iterations: float = 100,
            resolution: Literal["high", "medium",
                                "low"] = "medium") -> tasks.Task:
        """Simulates the wind tunnel scenario synchronously.

        Args:
            simulator: Simulator used to simulate the scenario.
                Valid simulators: OpenFOAM.
            object_path: Path to object inserted in the wind tunnel.
            num_iterations: Number of iterations of the simulator.
            resolution: Level of detail of the mesh used for the simulation.
                Options: "high", "medium" or "low".
            machine_group: The machine group to use for the simulation.
            storage_dir: The parent directory where simulation
            results will be stored.
        """
        simulator.override_api_method_prefix("wind_tunnel")

        if object_path is None:
            raise ValueError("WindTunnel is empty. Object path not specified.")

        self.object_path = utils.files.resolve_path(object_path)
        self.params["num_iterations"] = num_iterations
        self.params["resolution"] = MeshResolution[resolution.upper()].value

        commands = self.get_commands()
        task = super().simulate(simulator,
                                machine_group=machine_group,
                                storage_dir=storage_dir,
                                commands=commands)

        return task

    def add_extra_input_files(self, simulator: simulators.OpenFOAM, input_dir):  # pylint: disable=unused-argument
        """Configure object to be in specific place on input directory."""

        object_dir = os.path.join(input_dir, "constant", "triSurface")
        os.mkdir(object_dir)
        shutil.copy(self.object_path, os.path.join(object_dir, "object.obj"))
