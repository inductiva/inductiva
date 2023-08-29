"""Fluid tank scenario."""
from dataclasses import dataclass, field
from enum import Enum
from functools import singledispatchmethod
import os
from typing import List, Literal, Optional
from uuid import UUID

from inductiva import tasks, resources
from inductiva.scenarios import Scenario
from inductiva.simulation import Simulator
from inductiva.fluids.shapes import BaseShape
from inductiva.fluids.shapes import Rectangle
from inductiva.fluids.shapes import Circle
from inductiva.fluids.shapes import Cube
from inductiva.fluids.shapes import Cylinder
from inductiva.fluids.fluid_types import FluidType
from inductiva.fluids.fluid_types import WATER
from inductiva.fluids.simulators import SPlisHSPlasH
from inductiva.utils.templates import replace_params_in_template
from inductiva.fluids.scenarios.fluid_tank.output import FluidTankOutput

from . import mesh_file_utils

SPLISHSPLASH_TEMPLATE_FILENAME = "fluid_tank_template.splishsplash.json.jinja"
SPLISHSPLASH_CONFIG_FILENAME = "fluid_tank.json"
TANK_MESH_FILENAME = "tank.obj"
FLUID_MESH_FILENAME = "fluid.obj"


@dataclass
class ParticleRadius(Enum):
    """Sets particle radius according to resolution."""
    HIGH = 0.005
    MEDIUM = 0.01
    LOW = 0.02


@dataclass
class TimeStep(Enum):
    """Sets time step according to resolution."""
    HIGH = 0.0025
    MEDIUM = 0.005
    LOW = 0.01


# Tank inlets.
@dataclass
class BaseTankInlet:
    """Base tank inlet."""
    fluid_velocity: float = 1
    position: List[float] = field(default_factory=lambda: [0, 0])


@dataclass
class RectangularTankInlet(BaseTankInlet, Rectangle):
    """Rectangular tank inlet."""
    pass


@dataclass
class CircularTankInlet(BaseTankInlet, Circle):
    """Circular tank inlet."""
    pass


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

        self.shape = Cube(
            dimensions=dimensions,
            position=[*top_base_position, -dimensions[2]],
        )


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

        self.shape = Cylinder(
            radius=radius,
            height=height,
            position=[*top_base_position, -height],
        )


class FluidTank(Scenario):
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

    valid_simulators = [SPlisHSPlasH]

    def __init__(
        self,
        shape: BaseShape = Cylinder(radius=0.5, height=1),
        fluid: FluidType = WATER,
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
        simulator: Simulator = SPlisHSPlasH(),
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
        device: Literal["cpu", "gpu"] = "cpu",
        simulation_time: float = 5,
        resolution: Literal["low", "medium", "high"] = "low",
        output_time_step: float = 0.1,
        particle_sorting: bool = False,
    ) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: Simulator to use. Supported simulators are: SPlisHSPlasH.
            machine_group: The MachineGroup to use for the simulation.
            simulation_time: Total simulation time, in seconds.
            output_time_step: Time step for the output, in seconds.
            resolution: Resolution of the simulation. Controls the particle
              radius and time step. Accepted values are: "low", "medium",
              "high".
            particle_sorting: Whether to use particle sorting.
            run_async: Whether to run the simulation asynchronously.
        """

        self.simulation_time = simulation_time
        self.particle_radius = ParticleRadius[resolution.upper()].value
        self.time_step = TimeStep[resolution.upper()].value
        self.output_time_step = output_time_step
        self.particle_sorting = particle_sorting

        task = super().simulate(
            simulator,
            machine_group=machine_group,
            run_async=run_async,
            device=device,
            sim_config_filename=self.get_config_filename(simulator),
        )

        task.set_output_class(FluidTankOutput)

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

    def create_mesh_files(self, output_dir: str):
        """Creates mesh files for the tank and fluid in the given directory."""
        mesh_file_utils.create_tank_mesh_file(
            shape=self.shape,
            outlet=self.outlet,
            path=os.path.join(output_dir, TANK_MESH_FILENAME),
        )

        mesh_file_utils.create_tank_fluid_mesh_file(
            shape=self.shape,
            fluid_level=self.fluid_level,
            margin=2 * self.particle_radius,
            path=os.path.join(output_dir, FLUID_MESH_FILENAME),
        )

    @singledispatchmethod
    def get_config_filename(self, simulator: Simulator):
        pass

    @singledispatchmethod
    def create_input_files(self, simulator: Simulator):
        pass


@FluidTank.get_config_filename.register
def _(cls, simulator: SPlisHSPlasH) -> str:  # pylint: disable=unused-argument
    """Returns the config filename for SPlisHSPlasH."""
    return SPLISHSPLASH_CONFIG_FILENAME


@FluidTank.create_input_files.register
def _(self, simulator: SPlisHSPlasH, input_dir):  # pylint: disable=unused-argument
    """Creates SPlisHSPlasH simulation input files."""

    self.create_mesh_files(input_dir)

    bounding_box_min, bounding_box_max = self.get_bounding_box()
    inlet_position = [
        self.inlet.position[0],
        self.inlet.position[1],
        bounding_box_max[2],
    ]

    replace_params_in_template(
        template_path=os.path.join(os.path.dirname(__file__),
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
            "inlet_width": int(self.inlet.radius / self.particle_radius),
            "inlet_fluid_velocity": self.inlet.fluid_velocity,
            "bounding_box_min": bounding_box_min,
            "bounding_box_max": bounding_box_max,
        },
        output_file_path=os.path.join(input_dir, SPLISHSPLASH_CONFIG_FILENAME),
    )
