"""Wind terrain scenario for air flowing over complex terrains."""

from functools import singledispatchmethod
import io
import os
import shutil
from typing import List, Optional

import numpy as np

import inductiva
from inductiva import simulators, resources, scenarios, world, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "wind_terrain")
OPENFOAM_TEMPLATE_INPUT_DIR = "openfoam"
FILES_SUBDIR = "files"
COMMANDS_TEMPLATE_FILE_NAME = "commands.json.jinja"
TERRAIN_FILENAME = "terrain.stl"


class WindOverTerrain(scenarios.Scenario):
    """Wind flowing over complex terrain scenario.

    This simulation scenario models the steady-state conditions of
    wind flowing over complex terrain (e.g., mountains, valleys, etc.).

    The terrain is modeled through a 2D surface in a 3D world, defined
    through a mesh. The wind is injected through one of the side walls with
    a certain velocity vector and and at a certain region of the wall.
    Here, the wind is initialized only on a circular region of the wall
    to simulate gusts of wind.

    Schematic representation of the wind at the input wall and
    the simulation scenario both on a 2D xz-planexz (z is the vertical axis):

         Wind Injection                 2D Representation of Simulation
     ____________________                     ____________________
    |                    |                   |                    |
    |    * *             |                   |                    |
    |  *     *           |                   |->                  |
    |  *     *           |                   |->          /|      |
    |    * *             |                   |         /|/ |      |
    |____________________|                   |________/____|______|

    This scenario solves the steady-state continuity and momentum equations
    (time-independent) with the assumption of incompressible flow.
    The simulation solves the time-independent equations for several
    time steps, based on the state of the previous one. The end goal is
    to determine the steady-state of the system, i.e., where the flow
    does not change in time anymore.

    This scenario can be simulated with OpenFOAM.
    """

    valid_simulators = [simulators.OpenFOAM]

    @utils.optional_deps.needs_fluids_extra_deps
    def __init__(self,
                 terrain: world.Terrain,
                 wind_velocity: List[float],
                 wind_position: Optional[List[float]] = None,
                 atmosphere_height: float = 300):
        """Initializes the `WindTerrain` conditions.

        Args:
            wind_velocity (List): Velocity of the air flow (m/s).
            wind_position (List): Absolute Position of the wind flow (m).
                Note: The position needs to be above the terrain to occur
                any wind flow.
            terrain (inductiva.world.Terrain): Terrain object that describes
                the profile of the terrain with a mesh.
            atmosphere_height (float): Altitude (m) above the lowest point of
                terrain (m) that establishes the region of air where wind flows.
                Notice that the wind_position needs to be inside
                this atmosphere region.
        """

        self.wind_velocity = np.array(wind_velocity)
        # Compute the wind_direction from the velocity vector.
        self.wind_direction = self.wind_velocity / np.linalg.norm(
            self.wind_velocity)
        if wind_position is None:
            wind_position = [
                terrain.center["x"], terrain.center["y"],
                terrain.bounds["z"][1] + 10
            ]

        self.wind_position = wind_position
        self.terrain = terrain

        top_boundary_height = terrain.bounds["z"][0] + atmosphere_height

        if terrain.bounds["z"][1] >= top_boundary_height:
            raise ValueError(
                "The terrain height surpasses the top athmosphere boundary.")
        else:
            self.top_boundary_height = top_boundary_height

    def simulate(
        self,
        simulator: simulators.Simulator = simulators.OpenFOAM(),
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
        n_cores: int = 2,
        num_iterations: int = 100,
    ) -> inductiva.tasks.Task:
        """Simulates the wind over the terrain scenario.

        Args:
            simulator: Simulator used to simulate the scenario.
                Valid simulators: OpenFOAM.
            machine_group: The MachineGroup to use for the simulation.
            num_iterations: Number of iterations to run the simulation.
            n_cores: Number of cores to use for the simulation.
            run_async: Whether to run the simulation asynchronously.
        """
        simulator.override_api_method_prefix("wind_terrain")

        self.num_iterations = num_iterations
        self.n_cores = n_cores

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


@WindOverTerrain.create_input_files.register
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
            os.path.join("system", "blockMeshDict_template.openfoam.jinja"),
            os.path.join("system", "controlDict_template.openfoam.jinja"),
            os.path.join("system", "decomposeParDict_template.openfoam.jinja"),
            os.path.join("system", "snappyHexMeshDict_template.openfoam.jinja"),
            os.path.join("system", "topoSetDict_template.openfoam.jinja"),
            os.path.join("0", "include",
                         "ABLConditions_template.openfoam.jinja"),
            os.path.join("constant", "fvOptions_template.openfoam.jinja")
        ],
        params={
            "terrain_bounds": self.terrain.bounds,
            "terrain_center": self.terrain.center,
            "terrain_path": TERRAIN_FILENAME,
            "wind_flow_dir": self.wind_direction,
            "flow_position": self.wind_position,
            "top_boundary_height": self.top_boundary_height,
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
    terrain_file_path = os.path.join(terrain_dir, TERRAIN_FILENAME)
    self.terrain.to_text_file(terrain_file_path)
