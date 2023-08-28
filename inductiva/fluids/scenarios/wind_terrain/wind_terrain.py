"""Wind terrain scenario for air flowing over complex terrains."""

from functools import singledispatchmethod
import os
import shutil
import tempfile
from typing import List, Optional
import uuid

import numpy as np

import inductiva
from inductiva import fluids

SCENARIO_TEMPLATE_DIR = os.path.join(inductiva.utils.templates.TEMPLATES_PATH,
                                     "wind_terrain")
OPENFOAM_TEMPLATE_INPUT_DIR = "openfoam"
FILES_SUBDIR = "files"
COMMANDS_TEMPLATE_FILE_NAME = "commands.json.jinja"
TERRAIN_FILENAME = "terrain.stl"

class WindOverTerrain(inductiva.scenarios.Scenario):
    """Wind Flowing over complex terrain scenario.

    This simulation scenario models the steady-state conditions of
    wind flowing over complex terrain (e.g., mountains, valleys, etc.).
    
    The terrain is modeled through a 2D surface in a 3D world and defined
    through a mesh. The wind is injected through one of the side walls with
    a certain vectorail velocity and position. Here, the wind is initialized
    only on a circular region of the wall to simulate gusts of wind.

   Schematic representation of the Wind and Simulation scenario:

         Wind Injection                 2D Representation of Simulation
     ____________________                     ____________________
    |                    |                   |                    |
    |    * *             |                   |                    |
    |  *     *           |                   |->                  |
    |  *     *           |                   |->          /\      |
    |    * *             |                   |         /\/  \     |
    |____________________|                   |________/______\____|

    This scenario solves the steady-state continuity and momentum equations
    (time-independent) with the assumption of incompressible flow. 
    The simulation solves the time-independent equations for several
    time steps, based on the state of the previous one. The end goal is
    to determine the steady-state of the system, i.e., where the flow
    does not change in time anymore.

    This scenarion can be simulation with OpenFOAM.
    """

    valid_simulators = [fluids.simulators.OpenFOAM]

    def __init__(self,
                 terrain: inductiva.world.Terrain,
                 wind_velocity: List[float],
                 wind_position: Optional[List[float]] = None,
                 athmosphere_height: float = 300):
        """Initializes the `WindTerrain` conditions.

        Args:
            wind_velocity (List): Velocity of the air flow (m/s).
            wind_position (List): Position of the wind flow (m) above the height
            of the terrain.
            terrain (inductiva.world.Terrain): Terrain object that describes
                the profile of the terrain with a mesh.
            athmosphere_height (float): Altitude above the terrain (m) to establish an
                atmosphere. Notice that the wind_position needs to be inside the
                atmosphere region.
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

        self.top_boundary_height = terrain.bounds["z"][0] + athmosphere_height

        if terrain.bounds["z"][1] >= self.top_boundary_height:
            raise ValueError(
                "The terrain height surpasses the top athmosphere boundary.")

    def simulate(
        self,
        simulator: inductiva.simulation.Simulator = fluids.simulators.OpenFOAM(
        ),
        resource_pool_id: Optional[uuid.UUID] = None,
        run_async: bool = False,
        n_cores: int = 1,
        num_iterations: int = 100,
    ) -> inductiva.tasks.Task:
        """Simulates the wind tunnel scenario synchronously.

        Args:
            simulator: Simulator used to simulate the scenario.
                Valid simulators: OpenFOAM.
            num_iterations: Number of iterations to run the simulation.
            n_cores: Number of cores to use for the simulation.
            resource_pool_id: Id of the resource pool to use for the simulation.
                TODO: Change resource pool id to machine group.
            run_async: Whether to run the simulation asynchronously.
        """

        self.num_iterations = num_iterations
        self.n_cores = n_cores

        commands = self.get_commands()

        task = super().simulate(simulator,
                                resource_pool_id=resource_pool_id,
                                run_async=run_async,
                                n_cores=n_cores,
                                commands=commands)

        task.set_output_class(
            inductiva.fluids.post_processing.SteadyStateOutput)

        return task

    def get_commands(self):
        """Returns the commands for the simulation."""

        commands_template_path = os.path.join(SCENARIO_TEMPLATE_DIR,
                                              OPENFOAM_TEMPLATE_INPUT_DIR,
                                              COMMANDS_TEMPLATE_FILE_NAME)

        with tempfile.NamedTemporaryFile() as commands_file:
            inductiva.utils.templates.replace_params_in_template(
                template_path=commands_template_path,
                params={"n_cores": self.n_cores},
                output_file_path=commands_file.name,
            )

            commands = self.read_commands_from_file(commands_file.name)

        return commands

    @singledispatchmethod
    def create_input_files(self, simulator: inductiva.simulation.Simulator):
        pass


@WindTerrain.create_input_files.register
def _(self, simulator: fluids.simulators.OpenFOAM, input_dir):  # pylint: disable=unused-argument
    """Creates OpenFOAM simulation input files."""

    # The WindTunnel with OpenFOAM requires changing multiple files
    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      OPENFOAM_TEMPLATE_INPUT_DIR, FILES_SUBDIR)

    # Copy all files from the template dir to the input directory
    shutil.copytree(template_files_dir,
                    input_dir,
                    dirs_exist_ok=True,
                    symlinks=True)

    inductiva.utils.templates.batch_replace_params_in_template(
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
