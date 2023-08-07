"""WindTerrain scenario for air flowing over complex terrains."""

from dataclasses import dataclass
from enum import Enum
from functools import singledispatchmethod
import os
import shutil
import tempfile
from typing import Optional, List, Literal
from uuid import UUID

from absl import logging
import numpy as np
import pyvista as pv

from inductiva import tasks
from inductiva.types import Path
from inductiva.scenarios import Scenario
from inductiva.simulation import Simulator
from inductiva.fluids.simulators import OpenFOAM
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       batch_replace_params_in_template,
                                       replace_params_in_template)
from inductiva.utils import files
from inductiva.fluids.scenarios.wind_terrain import output


SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "wind_terrain")
OPENFOAM_TEMPLATE_INPUT_DIR = "openfoam"
FILES_SUBDIR = "files"
COMMANDS_TEMPLATE_FILE_NAME = "commands.json.jinja"


class Terrain:
    """Terrain object."""

    def __init__(self, terrain_file: str):
        self.path = terrain_file
        self.mesh = pv.read(terrain_file)

    @property
    def bounds(self):
        """Returns the bounds of the terrain."""

        terrain_bounds = {
            "x": [self.mesh.bounds[0], self.mesh.bounds[1]],
            "y": [self.mesh.bounds[2], self.mesh.bounds[3]],
            "z": [self.mesh.bounds[4], self.mesh.bounds[5]]
        }

        return terrain_bounds

    @property
    def center(self):
        """Returns the center of the terrain."""

        terrain_center = {
            "x": (self.mesh.bounds[0] + self.mesh.bounds[1]) / 2,
            "y": (self.mesh.bounds[2] + self.mesh.bounds[3]) / 2
        }

        return terrain_center


class WindTerrain(Scenario):
    """Physical scenario of a configurable wind terrain simulation.

    This scenario solves steady-state continuity and momentum equations
    (time-independent) with incompressible flow. 
    The simulation solves the time-independent equations for several
    time steps, based on the state of the previous one. The end goal is
    to determine the steady-state of the system, i.e., where the flow
    does not change in time anymore.
    """

    valid_simulators = [OpenFOAM]

    def __init__(self,
                 terrain: Terrain,
                 wind_velocity: List[float] = [10, 0, 0],
                 wind_position: List[float] = None,
                 altitude_m: float = 100):
        """Initializes the `WindTerrain` conditions.

        Args:
            wind_velocity (dict): Velocity of the air flow (m/s).
            terrain (Terrain): Terrain object with the profile of the terrain.
            altitude_m (float): Altitude above the terrain (m) to simulate
                the atmosphere.
        """


        self.wind_velocity = np.array(wind_velocity)
        self.wind_flow = self.wind_velocity / np.linalg.norm(self.wind_velocity)
        if wind_position is None:
            wind_position = [terrain.center["x"],
                             terrain.center["y"], 
                             terrain.bounds["z"][1] + 10]
            
        self.wind_position = wind_position
        self.terrain = terrain
        self.altitude_m = altitude_m

    def simulate(self,
                 simulator: Simulator = OpenFOAM(),
                 num_iterations: int = 100,
                 resource_pool_id: Optional[UUID] = None,
                 run_async: bool = False,
                 n_cores: int = 1) -> tasks.Task:
        """Simulates the wind tunnel scenario synchronously.

        Args:
            simulator: Simulator used to simulate the scenario.
                Valid simulators: OpenFOAM.
            num_iterations: Number of iterations to run the simulation.
            n_cores: Number of cores to use for the simulation.
            resource_pool_id: Id of the resource pool to use for the simulation.
        """

        self.num_iterations = num_iterations
        self.n_cores = n_cores

        commands = self.get_commands()

        task = super().simulate(simulator,
                                resource_pool_id=resource_pool_id,
                                run_async=run_async,
                                n_cores=n_cores,
                                commands=commands)

        return task

    def get_commands(self):
        """Returns the commands for the simulation."""

        commands_template_path = os.path.join(SCENARIO_TEMPLATE_DIR,
                                              OPENFOAM_TEMPLATE_INPUT_DIR,
                                              COMMANDS_TEMPLATE_FILE_NAME)

        with tempfile.NamedTemporaryFile() as commands_file:
            replace_params_in_template(
                template_path=commands_template_path,
                params={"n_cores": self.n_cores},
                output_file_path=commands_file.name,
            )

            commands = self.read_commands_from_file(commands_file.name)

        return commands

    @singledispatchmethod
    def create_input_files(self, simulator: Simulator):
        pass


@WindTerrain.create_input_files.register
def _(self, simulator: OpenFOAM, input_dir):  # pylint: disable=unused-argument
    """Creates OpenFOAM simulation input files."""

    # The WindTunnel with OpenFOAM requires changing multiple files
    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      OPENFOAM_TEMPLATE_INPUT_DIR, FILES_SUBDIR)

    # Copy all files from the template dir to the input directory
    shutil.copytree(template_files_dir,
                    input_dir,
                    dirs_exist_ok=True,
                    symlinks=True)

    batch_replace_params_in_template(
        templates_dir=input_dir,
        template_filenames=[
            os.path.join("system", "blockMeshDict_template.openfoam.jinja"),
            os.path.join("system", "controlDict_template.openfoam.jinja"),
            os.path.join("system", "decomposeParDict_template.openfoam.jinja"),
            os.path.join("system", "snappyHexMeshDict_template.openfoam.jinja"),
            os.path.join("system", "topoSetDict_template.openfoam.jinja"),
            os.path.join("0", "include", "ABLConditions_template.openfoam.jinja"),
            os.path.join("constant", "fvOptions_template.openfoam.jinja")
        ],
        params={
            "terrain_bounds": self.terrain.bounds,
            "terrain_center": self.terrain.center,
            "terrain_path": self.terrain.path,
            "wind_flow_dir": self.wind_flow,
            "flow_position": self.wind_position,
            "altitude": self.altitude_m,
            "num_iterations": self.num_iterations,
            "n_cores": self.n_cores,
        },
        output_filename_paths=[
            os.path.join(input_dir, "system", "blockMeshDict"),
            os.path.join(input_dir, "system", "controlDict"),
            os.path.join(input_dir, "system", "decomposeParDict"),
            os.path.join(input_dir, "system", "snappyHexMeshDict"),
            os.path.join(input_dir, "system", "topoSetDict"),
            os.path.join(input_dir, "0", "include", "ABLConditions"),
            os.path.join(input_dir, "constant", "fvOptions")
        ],
        remove_templates=True,
    )

    # Add terrain path to its respective place
    terrain_dir = os.path.join(input_dir, "constant", "triSurface")
    os.mkdir(terrain_dir)
    shutil.copy(self.terrain.path, os.path.join(terrain_dir, self.terrain.path))
