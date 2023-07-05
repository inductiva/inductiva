"""Classes describing fluid tank scenarios."""

from dataclasses import dataclass, field
from enum import Enum
from functools import singledispatchmethod
import os
from typing import List, Literal, Optional
from uuid import UUID

from inductiva.types import Path
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
from inductiva.fluids.post_processing.splishsplash import convert_vtk_data_dir_to_netcdf
from inductiva.utils.templates import replace_params_in_template

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
    """Fluid tank."""

    def __init__(
        self,
        shape: BaseShape = Cylinder(radius=0.5, height=1),
        fluid: FluidType = WATER,
        fluid_level: float = 0,
        inlet: Optional[BaseTankInlet] = CircularTankInlet(radius=0.1),
        outlet: Optional[BaseTankOutlet] = CylindricalTankOutlet(radius=0.1,
                                                                 height=0.1),
    ):
        """Initializes a fluid tank.

        Args:
            shape: Shape of the tank.
            fluid: Fluid type.
            fluid_level: Fluid level initially in the tank.
            inlet: Inlet of the tank.
            outlet: Outlet of the tank.
        """
        self.shape = shape
        self.fluid = fluid
        self.fluid_level = fluid_level
        self.inlet = inlet
        self.outlet = outlet

    def simulate(
        self,
        simulator: Simulator = SPlisHSPlasH(),
        output_dir: Optional[Path] = None,
        resource_pool_id: Optional[UUID] = None,
        device: Literal["cpu", "gpu"] = "cpu",
        simulation_time: float = 5,
        resolution: Literal["low", "medium", "high"] = "low",
        output_time_step: float = 0.1,
        particle_sorting: bool = False,
    ):
        """Simulates the scenario.

        Args:
            simulator: Simulator to use.
            output_dir: Directory to store the simulation output.
            device: Device in which to run the simulation.
            simulation_time: Simulation time, in seconds.
            output_time_step: Time step between output files, in seconds.
            resolution: Resolution of the simulation. Controls the particle
                radius and time step.
            particle_sorting: Whether to use particle sorting.
        """

        self.simulation_time = simulation_time
        self.particle_radius = ParticleRadius[resolution.upper()].value
        self.time_step = TimeStep[resolution.upper()].value
        self.output_time_step = output_time_step
        self.particle_sorting = particle_sorting

        output_path = super().simulate(
            simulator,
            output_dir=output_dir,
            resource_pool_id=resource_pool_id,
            device=device,
        )

        # TODO: Add any kind of post-processing here, e.g. convert files?
        convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(output_path, "vtk"),
            output_time_step=self.output_time_step,
            netcdf_data_dir=os.path.join(output_path, "netcdf"))

        return output_path

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


@FluidTank.get_config_filename.register
def _(cls, simulator: SPlisHSPlasH) -> str:  # pylint: disable=unused-argument
    """Returns the config filename for SPlisHSPlasH."""
    return SPLISHSPLASH_CONFIG_FILENAME


@FluidTank.gen_aux_files.register
def _(self, simulator: SPlisHSPlasH, input_dir):  # pylint: disable=unused-argument
    """Generates auxiliary files for SPlisHSPlasH."""
    mesh_file_utils.create_tank_mesh_file(
        shape=self.shape,
        outlet=self.outlet,
        path=os.path.join(input_dir, TANK_MESH_FILENAME),
    )

    mesh_file_utils.create_tank_fluid_mesh_file(
        shape=self.shape,
        fluid_level=self.fluid_level,
        margin=2 * self.particle_radius,
        path=os.path.join(input_dir, FLUID_MESH_FILENAME),
    )


@FluidTank.gen_config.register
def _(self, simulator: SPlisHSPlasH, input_dir: str):  # pylint: disable=unused-argument
    """Generates the configuration file for SPlisHSPlasH."""

    bounding_box_min, bounding_box_max = self.get_bounding_box()
    inlet_position = [
        self.inlet.position[0],
        self.inlet.position[1],
        bounding_box_max[2],
    ]

    replace_params_in_template(
        templates_dir=os.path.dirname(__file__),
        template_filename=SPLISHSPLASH_TEMPLATE_FILENAME,
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
