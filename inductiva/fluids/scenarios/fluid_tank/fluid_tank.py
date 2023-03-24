"""Classes describing fluid tank scenarios."""

from dataclasses import dataclass, field
import os
import tempfile
from typing import List, Literal, Optional, Union

from inductiva.types import Path
from inductiva.fluids.shapes import BaseShape
from inductiva.fluids.shapes import Rectangle
from inductiva.fluids.shapes import Circle
from inductiva.fluids.shapes import Cube
from inductiva.fluids.shapes import Cylinder
from inductiva.fluids.fluid_types import FluidType
from inductiva.fluids.fluid_types import WATER
from inductiva.fluids.simulators import SPlisHSPlasH
from inductiva.fluids.simulators import SPlisHSPlasHParameters
from inductiva.fluids.simulators import DualSPHysicsParameters
from inductiva.fluids.post_processing.splishsplash import convert_vtk_data_dir_to_netcdf
from inductiva.utils.templates import replace_params_in_template

from inductiva.fluids._output_post_processing import SimulationOutput

from . import mesh_file_utils

SPLISHSPLASH_TEMPLATE_FILENAME = "fluid_tank_template.splishsplash.json.jinja"
SPLISHSPLASH_INPUT_FILENAME = "fluid_tank.json"
TANK_MESH_FILENAME = "tank.obj"
FLUID_MESH_FILENAME = "fluid.obj"

SIMULATION_TIME = 5
OUTPUT_TIME_STEP = 0.1
TIME_STEP = 0.005
PARTICLE_RADIUS = 0.02
Z_SORT = True


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


class FluidTank:
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
        device: Literal["cpu", "gpu"] = "cpu",
        engine: Literal["SPlisHSPlasH"] = "SPlisHSPlasH",
        engine_params: Union[DualSPHysicsParameters,
                             SPlisHSPlasHParameters] = SPlisHSPlasHParameters(),
        output_dir: Optional[Path] = None,
    ):
        """Simulates the fluid tank.
        
        Args:
            device: Device to run the simulation on.
            engine: Simulation engine to use.
            engine_params: Engine parameters.
            output_dir: Directory to store simulation output.

        Returns:
            Simulation output.
        """

        # Create a temporary directory to store simulation input files
        self.input_temp_dir = tempfile.TemporaryDirectory()  #pylint: disable=consider-using-with

        if engine.lower() == "splishsplash":
            sim_output_path = self._simulate_with_splishsplash(
                device=device,
                engine_params=engine_params,
                output_dir=output_dir)
        elif engine.lower() == "dualsphysics":
            raise NotImplementedError(
                "The engine 'DualSPHysics' is not supported yet "
                "for fluid tank simulations.")
        else:
            raise ValueError(f"Invalid engine `{engine}`.")

        # Delete temporary input directory
        self.input_temp_dir.cleanup()

        return SimulationOutput(sim_output_path)

    def _simulate_with_splishsplash(self, device, engine_params, output_dir):
        """Simulates the fluid tank with SPlisHSPlasH."""

        if not isinstance(engine_params, SPlisHSPlasHParameters):
            raise ValueError(f"Invalid engine parameters `{engine_params}`.")

        input_dir = self.input_temp_dir.name

        self._create_splishsplash_aux_files(input_dir)
        self._replace_params_in_splishsplash_template(input_dir)

        simulator = SPlisHSPlasH(sim_dir=input_dir,
                                 input_filename=SPLISHSPLASH_INPUT_FILENAME)

        output_path = simulator.simulate(device=device, output_dir=output_dir)

        convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(output_path, "vtk"),
            output_time_step=OUTPUT_TIME_STEP,
            netcdf_data_dir=os.path.join(output_path, "netcdf"))

        return output_path

    def _create_splishsplash_aux_files(self, input_dir):
        """Creates auxiliary files for SPlisHSPlasH simulation."""

        mesh_file_utils.create_tank_mesh_file(
            shape=self.shape,
            outlet=self.outlet,
            path=os.path.join(input_dir, TANK_MESH_FILENAME),
        )

        mesh_file_utils.create_tank_fluid_mesh_file(
            shape=self.shape,
            fluid_level=self.fluid_level,
            margin=2 * PARTICLE_RADIUS,
            path=os.path.join(input_dir, FLUID_MESH_FILENAME),
        )

    def _replace_params_in_splishsplash_template(self, input_dir):
        """Replaces parameters in SPlisHSPlasH template input file."""

        bounding_box_min, bounding_box_max = self._get_bounding_box()
        inlet_position = [
            self.inlet.position[0],
            self.inlet.position[1],
            bounding_box_max[2],
        ]

        replace_params_in_template(
            templates_dir=os.path.dirname(__file__),
            template_filename=SPLISHSPLASH_TEMPLATE_FILENAME,
            params={
                "simulation_time": SIMULATION_TIME,
                "time_step": TIME_STEP,
                "particle_radius": PARTICLE_RADIUS,
                "data_export_rate": 1 / OUTPUT_TIME_STEP,
                "z_sort": Z_SORT,
                "tank_filename": TANK_MESH_FILENAME,
                "fluid_filename": FLUID_MESH_FILENAME,
                "fluid": self.fluid,
                "inlet_position": inlet_position,
                "inlet_width": int(self.inlet.radius / PARTICLE_RADIUS),
                "inlet_fluid_velocity": self.inlet.fluid_velocity,
                "bounding_box_min": bounding_box_min,
                "bounding_box_max": bounding_box_max,
            },
            output_file_path=os.path.join(input_dir,
                                          SPLISHSPLASH_INPUT_FILENAME),
        )

    def _get_bounding_box(self):
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
