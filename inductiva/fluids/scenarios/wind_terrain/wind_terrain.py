"""Wind terrain scenario for air flowing over complex terrains."""

from functools import singledispatchmethod
import os
import random
import shutil
import tempfile
import typing
import uuid

import numpy as np
import pyvista as pv

import inductiva
from inductiva import fluids

SCENARIO_TEMPLATE_DIR = os.path.join(inductiva.utils.templates.TEMPLATES_PATH,
                                     "wind_terrain")
OPENFOAM_TEMPLATE_INPUT_DIR = "openfoam"
FILES_SUBDIR = "files"
COMMANDS_TEMPLATE_FILE_NAME = "commands.json.jinja"
TERRAIN_FILENAME = "terrain.stl"


class Terrain:
    """Terrain object."""

    def __init__(self, terrain_mesh):
        """Initializes the `Terrain` object.
        
        Args:
            terrain_mesh (pyvista.StructuredGrid): Terrain mesh."""
        self.mesh = terrain_mesh

    @classmethod
    def from_file(cls, terrain_file: str):
        """Creates a `Terrain` object from a text file."""
        terrain_mesh = pv.read(terrain_file)

        return cls(terrain_mesh)

    @classmethod
    def from_random_generation(cls,
                               x_range: typing.Sequence[float],
                               y_range: typing.Sequence[float],
                               x_num: int,
                               y_num: int,
                               height_factor: float = 10,
                               initial_roughness: float = 1,
                               roughness_factor: float = 0.5):
        """Creates a `Terrain` object with random elevations.
        
        The elevation of the corners are chosen randomly in the 
        range from 0 to 1. This limits the maximum height of the terrain
        due to the nature of the terrain generation algorithm.
        To increase the maximum height we multiply by the `height_factor`
        to obtain the terrain

        TODO: Improve terrain generation for high discrepancy
        terrains.
        """
        corner_values = [
            random.uniform(0, 1),
            random.uniform(0, 1),
            random.uniform(0, 1),
            random.uniform(0, 1)
        ]

        x_grid, y_grid, z_elevation = inductiva.generative.generate_grid_elevation(
            x_range, y_range, x_num, y_num, corner_values, initial_roughness,
            roughness_factor)

        terrain = pv.StructuredGrid(x_grid, y_grid, z_elevation * height_factor)

        return cls(terrain)

    def to_text_file(self, text_file_path: str):
        """Saves the terrain to a text file.
        
        Args
            text_file_path: Path to the text file
        """

        terrain_geometry = self.mesh.extract_geometry()
        terrain_geometry.save(text_file_path)

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

    @property
    def mesh_resolution(self):
        """Returns the resolution of the terrain mesh."""

        num_cells = self.mesh.n_cells
        num_points = self.mesh.n_points

        return num_cells, num_points

    def plot(self):
        """Returns a plot of the terrain."""

        plot = pv.Plotter()
        plot.add_mesh(self.mesh)
        plot.show()
        plot.close()


class WindTerrain(inductiva.scenarios.Scenario):
    """Physical scenario of a configurable wind terrain simulation.

    This scenario solves steady-state continuity and momentum equations
    (time-independent) with incompressible flow. 
    The simulation solves the time-independent equations for several
    time steps, based on the state of the previous one. The end goal is
    to determine the steady-state of the system, i.e., where the flow
    does not change in time anymore.
    """

    valid_simulators = [fluids.simulators.OpenFOAM]

    def __init__(self,
                 terrain: Terrain,
                 wind_velocity: typing.List[float],
                 wind_position: typing.Optional[typing.List[float]] = None,
                 altitude_m: float = 300):
        """Initializes the `WindTerrain` conditions.

        Args:
            wind_velocity (List): Velocity of the air flow (m/s).
            wind_position (List): Position of the wind flow (m) above the height
            of the terrain.
            terrain (Terrain): Terrain object with the profile of the terrain.
            altitude_m (float): Altitude above the terrain (m) to establish an
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
        self.altitude_m = altitude_m

    def simulate(
        self,
        simulator: inductiva.simulation.Simulator = fluids.simulators.OpenFOAM(),
        resource_pool_id: typing.Optional[uuid.UUID] = None,
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
        """

        self.num_iterations = num_iterations
        self.n_cores = n_cores

        commands = self.get_commands()

        task = super().simulate(simulator,
                                resource_pool_id=resource_pool_id,
                                run_async=run_async,
                                n_cores=n_cores,
                                commands=commands)

        task.set_output_class(inductiva.fluids.post_processing.SteadyStateOutput)

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
    terrain_file_path = os.path.join(terrain_dir, TERRAIN_FILENAME)
    self.terrain.to_text_file(terrain_file_path)
