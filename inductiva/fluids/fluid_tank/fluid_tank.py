"""Fluid tank scenario."""
from dataclasses import dataclass
import enum
from functools import singledispatchmethod
import json
import os
from typing import List, Literal, Optional

from inductiva import resources, fluids, simulators, scenarios, tasks, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "fluid_tank")
SPLISHSPLASH_TEMPLATE_INPUT_DIR = "splishsplash"
SPLISHSPLASH_TEMPLATE_FILENAME = "fluid_tank_template.splishsplash.json.jinja"
SPLISHSPLASH_CONFIG_FILENAME = "fluid_tank.json"
TANK_JSON_FILENAME = "tank.json"
TANK_MESH_FILENAME = "tank.obj"
FLUID_MESH_FILENAME = "fluid.obj"


@dataclass
class ParticleRadius(enum.Enum):
    """Sets particle radius according to resolution."""
    HIGH = 0.005
    MEDIUM = 0.01
    LOW = 0.02


@dataclass
class TimeStep(enum.Enum):
    """Sets time step according to resolution."""
    HIGH = 0.0025
    MEDIUM = 0.005
    LOW = 0.01


# Tank inlets.
class BaseTankInlet:
    """Base tank inlet."""

    def __init__(self, fluid_velocity: float):
        """Initializes a base tank inlet.

        Args:
            fluid_velocity: Fluid velocity.
        """
        self.fluid_velocity = fluid_velocity

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the inlet."""
        return {
            "fluid_velocity": self.fluid_velocity,
        }


class CircularTankInlet(BaseTankInlet):
    """Circular tank inlet."""

    def __init__(self,
                 fluid_velocity: float = 1,
                 position: List[float] = (0, 0),
                 radius: float = 0.1):
        """Initializes a circular tank inlet."""
        self.shape = fluids.shapes.Circle(radius=radius, position=position)
        super().__init__(fluid_velocity=fluid_velocity)

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the inlet."""
        return {
            **super().to_dict(),
            "shape": self.shape.to_dict(),
        }


# Tank outlets.
class BaseTankOutlet:
    """Base tank outlet."""
    pass


class CubicTankOutlet(BaseTankOutlet):
    """Cubic tank outlet."""

    def __init__(
        self,
        dimensions: List[float],
        top_base_position: Optional[List[float]] = None,
    ):
        """Initializes a cubic tank outlet.

        Args:
            dimensions: Dimensions of the outlet.
            top_base_position: Position of the top base of the outlet.
        """
        top_base_position = top_base_position or [0, 0]

        self.shape = fluids.shapes.Cube(
            dimensions=dimensions,
            position=[*top_base_position, -dimensions[2]],
        )

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the outlet."""
        return {"shape": self.shape.to_dict()}


class CylindricalTankOutlet(BaseTankOutlet):
    """Cylindrical tank outlet."""

    def __init__(
        self,
        radius: float,
        height: float,
        top_base_position: Optional[List[float]] = None,
    ):
        """Initializes a cylindrical tank outlet.

        Args:
            radius: Radius of the outlet.
            height: Height of the outlet.
            top_base_position: Position of the top base of the outlet.
        """
        top_base_position = top_base_position or [0, 0]

        self.shape = fluids.shapes.Cylinder(
            radius=radius,
            height=height,
            position=[*top_base_position, -height],
        )

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the outlet."""
        return {"shape": self.shape.to_dict()}


class FluidTank(scenarios.Scenario):
    """Fluid tank scenario.

    This is a simulation scenario for a fluid tank. The tank has a 3D shape that
    may be cubic or cylindrical. Fluid is injected in the tank via an inlet
    located at the top of the tank, and flows out of the tank via an outlet
    located at the bottom of the tank. The motion of the fluid is controled
    by gravity. The fluid properties such as density and kinematic viscosity
    are configurable. The initial fluid level in the tank is also configurable,
    as well as the inlet and outlet positions and dimensions.

    The main axis of the tank is the z axis. The inlet and outlet are located at
    the top and bottom of the tank, respectively, along the z axis. Fluid is
    injected from the inlet in the negative z direction with a given velocity.

    Schematic representation of the simulation scenario: e.g. x/y points right,
    z points up.

       inlet
    _____________________
    |    |              |
    |    v              |
    |                   |
    |                   |
    |___________________|  fluid level
    |                   |
    |                   |
    |                   |
    |                   |
    |___________________|
                    |
                    v  outlet

    The scenario can be simulated with SPlisHSPlasH.
    """

    valid_simulators = [simulators.SplishSplash]

    def __init__(
        self,
        shape: fluids.shapes.BaseShape = fluids.shapes.Cylinder(radius=0.5,
                                                                height=1),
        fluid: fluids.FluidType = fluids.WATER,
        fluid_level: float = 0.5,
        inlet: Optional[BaseTankInlet] = CircularTankInlet(radius=0.1),
        outlet: Optional[BaseTankOutlet] = CylindricalTankOutlet(radius=0.1,
                                                                 height=0.1),
    ):
        """Initializes a fluid tank.

        Args:
            shape: The shape of the tank.
            fluid: The fluid type.
            fluid_level: The fluid level initially in the tank.
            inlet: The inlet of the tank.
            outlet: The outlet of the tank.
        """
        self.shape = shape
        self.fluid = fluid
        self.fluid_level = fluid_level
        self.inlet = inlet
        self.outlet = outlet

    def simulate(
        self,
        simulator: simulators.Simulator = simulators.SplishSplash(),
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
        simulation_time: float = 1,
        resolution: Literal["low", "medium", "high"] = "low",
        output_time_step: float = 0.1,
        particle_sorting: bool = False,
    ) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: Simulator to use. Supported simulators are: SPlisHSPlasH.
            machine_group: The machine group to use for the simulation.
            simulation_time: Total simulation time, in seconds.
            output_time_step: Time step for the output, in seconds.
            resolution: Resolution of the simulation. Controls the particle
              radius and time step. Accepted values are: "low", "medium",
              "high".
            particle_sorting: Whether to use particle sorting.
            run_async: Whether to run the simulation asynchronously.
        """
        simulator.override_api_method_prefix("fluid_tank")

        self.simulation_time = simulation_time
        self.particle_radius = ParticleRadius[resolution.upper()].value
        self.time_step = TimeStep[resolution.upper()].value
        self.output_time_step = output_time_step
        self.particle_sorting = particle_sorting

        task = super().simulate(
            simulator,
            machine_group=machine_group,
            run_async=run_async,
            particle_radius=self.particle_radius,
            sim_config_filename=self.get_config_filename(simulator),
        )

        return task

    def get_bounding_box(self):
        """Gets the bounding box of the tank.

        Returns:
            Tuple of two lists representing the minimum and maximum coordinates
            of the bounding box of the tank, respectively.
        """

        bounding_box_min, bounding_box_max = self.shape.get_bounding_box()

        # Extend the bounding box to include the outlet.
        if self.outlet is not None:
            outlet_bounding_box_min, _ = self.outlet.shape.get_bounding_box()
            bounding_box_min[2] = outlet_bounding_box_min[2]

        return bounding_box_min, bounding_box_max

    def to_dict(self) -> dict:
        """Returns a dictionary representation of the scenario."""
        return {
            "shape": self.shape.to_dict(),
            "fluid": self.fluid.to_dict(),
            "fluid_level": self.fluid_level,
            "inlet": self.inlet.to_dict(),
            "outlet": self.outlet.to_dict(),
        }

    def create_json_file(self, output_path):
        """Creates a JSON file with the scenario parameters."""

        with open(output_path, "w", encoding="utf-8") as file:
            json.dump(self.to_dict(), file)

    @singledispatchmethod
    def get_config_filename(self, simulator: simulators.Simulator):
        pass

    @singledispatchmethod
    def create_input_files(self, simulator: simulators.Simulator):
        pass


@FluidTank.get_config_filename.register
def _(cls, simulator: simulators.SplishSplash) -> str:  # pylint: disable=unused-argument
    """Returns the config filename for SPlisHSPlasH."""
    return SPLISHSPLASH_CONFIG_FILENAME


@FluidTank.create_input_files.register
def _(self, simulator: simulators.SplishSplash, input_dir):  # pylint: disable=unused-argument
    """Creates SPlisHSPlasH simulation input files."""

    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      SPLISHSPLASH_TEMPLATE_INPUT_DIR)
    self.create_json_file(os.path.join(input_dir, TANK_JSON_FILENAME))

    bounding_box_min, bounding_box_max = self.get_bounding_box()
    inlet_position = [
        self.inlet.shape.position[0],
        self.inlet.shape.position[1],
        bounding_box_max[2],
    ]

    utils.templates.replace_params(
        template_path=os.path.join(template_files_dir,
                                   SPLISHSPLASH_TEMPLATE_FILENAME),
        params={
            "simulation_time": self.simulation_time,
            "time_step": self.time_step,
            "particle_radius": self.particle_radius,
            "data_export_rate": 1 / self.output_time_step,
            "z_sort": self.particle_sorting,
            "tank_filename": TANK_MESH_FILENAME,
            "fluid_filename": FLUID_MESH_FILENAME,
            "fluid": self.fluid,
            "inlet_position": inlet_position,
            "inlet_width": int(self.inlet.shape.radius / self.particle_radius),
            "inlet_fluid_velocity": self.inlet.fluid_velocity,
            "bounding_box_min": bounding_box_min,
            "bounding_box_max": bounding_box_max,
        },
        output_file=os.path.join(input_dir, SPLISHSPLASH_CONFIG_FILENAME),
    )
