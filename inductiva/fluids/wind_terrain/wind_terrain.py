"""Wind terrain scenario for air flowing over complex terrains."""
import os
from typing import List, Optional

import numpy as np

import inductiva
from inductiva import simulators, resources, scenarios, world, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "wind_terrain")
OPENFOAM_TEMPLATE_INPUT_DIR = "openfoam"


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
    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      OPENFOAM_TEMPLATE_INPUT_DIR)

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
        wind_velocity = np.array(wind_velocity)
        self.params["wind_velocity"] = wind_velocity
        # Compute the wind direction from the velocity vector.
        self.params["wind_flow_dir"] = wind_velocity / np.linalg.norm(
            wind_velocity)
        if wind_position is None:
            wind_position = [
                terrain.center["x"], terrain.center["y"],
                terrain.bounds["z"][1] + 10
            ]

        self.params["flow_position"] = wind_position
        self.params["terrain"] = terrain

        top_boundary_height = terrain.bounds["z"][0] + atmosphere_height

        if terrain.bounds["z"][1] >= top_boundary_height:
            raise ValueError(
                "The terrain height surpasses the top athmosphere boundary.")

        self.params["top_boundary_height"] = top_boundary_height

    def simulate(
        self,
        simulator: simulators.Simulator = simulators.OpenFOAM(),
        machine_group: Optional[resources.MachineGroup] = None,
        storage_dir: Optional[str] = "",
        num_iterations: int = 100,
    ) -> inductiva.tasks.Task:
        """Simulates the wind over the terrain scenario.

        Args:
            simulator: Simulator used to simulate the scenario.
                Valid simulators: OpenFOAM.
            machine_group: The MachineGroup to use for the simulation.
            num_iterations: Number of iterations to run the simulation.
            storage_dir: The parent directory where simulation
            results will be stored.
        """
        simulator.override_api_method_prefix("wind_terrain")

        self.params["num_iterations"] = num_iterations

        commands = self.get_commands()

        task = super().simulate(simulator,
                                machine_group=machine_group,
                                storage_dir=storage_dir,
                                commands=commands)

        return task

    def add_extra_input_files(self, simulator: simulators.OpenFOAM, input_dir):  # pylint: disable=unused-argument
        """Configure object to be in specific place on input directory."""

        terrain_dir = os.path.join(input_dir, "constant", "triSurface")
        os.mkdir(terrain_dir)
        terrain_file_path = os.path.join(terrain_dir, "terrain.stl")
        self.params["terrain"].to_text_file(terrain_file_path)
