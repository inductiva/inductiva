"""Classes describing fluid tank scenarios."""

from dataclasses import dataclass, field
import os
import tempfile
from typing import List, Literal, Optional, Union

from inductiva_sph.splishsplash.io_utils import convert_vtk_data_dir_to_netcdf

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
from inductiva.utils.templates import replace_params_in_template_file

from inductiva.fluids._output_post_processing import SimulationOutput

from . import gmsh_utils

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
@dataclass
class BaseTankOutlet:
    """Base tank outlet."""
    top_base_position: List[float] = field(default_factory=lambda: [0, 0])


@dataclass
class CubicTankOutlet(BaseTankOutlet, Cube):
    """Cubic tank outlet."""
    pass


@dataclass
class CylindricalTankOutlet(BaseTankOutlet, Cylinder):
    """Cylindrical tank outlet."""
    pass


class FluidTank:
    """Fluid tank."""

    def __init__(
        self,
        shape: BaseShape = Cube(dimensions=[1, 1, 1]),
        fluid: FluidType = WATER,
        fluid_level: float = 0,
        inlet: Optional[BaseTankInlet] = CircularTankInlet(
            radius=0.1,
            position=[0, 0],
        ),
        outlet: Optional[BaseTankOutlet] = CylindricalTankOutlet(
            radius=0.1,
            height=0.1,
            top_base_position=[0, 0],
        ),
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

        _create_tank_mesh_file(
            shape=self.shape,
            outlet=self.outlet,
            path=os.path.join(input_dir, TANK_MESH_FILENAME),
        )

        _create_fluid_mesh_file(
            shape=self.shape,
            fluid_level=self.fluid_level,
            margin=2 * PARTICLE_RADIUS,
            path=os.path.join(input_dir, FLUID_MESH_FILENAME),
        )

        bounding_box_min, bounding_box_max = self._get_bounding_box()
        inlet_position = [
            self.inlet.position[0],
            self.inlet.position[1],
            bounding_box_max[2],
        ]

        replace_params_in_template_file(
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

        simulator = SPlisHSPlasH(sim_dir=input_dir,
                                 input_filename=SPLISHSPLASH_INPUT_FILENAME)

        output_path = simulator.simulate(device=device, output_dir=output_dir)

        convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(output_path, "vtk"),
            output_time_step=OUTPUT_TIME_STEP,
            netcdf_data_dir=os.path.join(output_path, "netcdf"))

        return output_path

    def _get_bounding_box(self):
        """Gets the bounding box of the tank.
        
        Returns:
            Tuple of two lists representing the minimum and maximum coordinates
            of the bounding box of the tank, respectively.
        """

        # Get the bounding box of the tank.
        if isinstance(self.shape, Cylinder):
            bounding_box_min = [
                -self.shape.radius,
                -self.shape.radius,
                0,
            ]
            bounding_box_max = [
                self.shape.radius,
                self.shape.radius,
                self.shape.height,
            ]

        elif isinstance(self.shape, Cube):
            bounding_box_min = [
                -self.shape.dimensions[0] / 2,
                -self.shape.dimensions[1] / 2,
                0,
            ]
            bounding_box_max = [
                self.shape.dimensions[0] / 2,
                self.shape.dimensions[1] / 2,
                self.shape.dimensions[2],
            ]
        else:
            raise ValueError(f"Invalid tank shape `{self.shape}`.")

        # Extend the bounding box to include the outlet.
        if self.outlet is not None:
            if isinstance(self.outlet, Cylinder):
                bounding_box_min[2] = -self.outlet.height
            elif isinstance(self.outlet, Cube):
                bounding_box_min[2] = -self.outlet.dimensions[2]
            else:
                raise ValueError(f"Invalid outlet shape `{self.outlet}`.")

        return bounding_box_min, bounding_box_max


def _create_tank_mesh_file(shape, outlet, path: str):
    """Creates a mesh file for the tank.
    
    The tank is composed of two blocks:
    - a main (cylindrical/cubic) block representing the tank itself;
    - an optional smaller (cylindrical/cubic) block representing a fluid outlet.
      When present, the top base of this block connects with the bottom base of
      the tank, such that fluid flows freely from the tank to the outlet. The
      bottom base of the outlet is also open, such that flow exits the outlet.
    
    Both blocks are assumed to have their main axes aligned with the z
    axis.
    
    Args:
        shape: Shape of the tank.
        outlet: Shape of the outlet. If `None`, no outlet is present.
        path: Path of the file to be created.
    """

    with gmsh_utils.gmshAPIWrapper():
        tank_base_hole_loops = []

        if outlet is not None:

            # Add a circle arc/rectangle loop representing the top base of the
            # outlet. An arc/loop is used instead of a circle/rectangle because
            # this face is not filled, i.e. it is not a surface.
            if isinstance(outlet, Cylinder):
                p_top_outlet, c_top_outlet, l_top_outlet = \
                    gmsh_utils.add_circle_arc(
                        x=outlet.top_base_position[0],
                        y=outlet.top_base_position[1],
                        z=0,
                        r=outlet.radius,
                    )
            elif isinstance(outlet, Cube):
                p_top_outlet, c_top_outlet, l_top_outlet = \
                    gmsh_utils.add_z_rectangle_loop(
                        x=outlet.top_base_position[0],
                        y=outlet.top_base_position[1],
                        z=0,
                        lx=outlet.dimensions[0],
                        ly=outlet.dimensions[1],
                    )

            # Add a circle arc/rectangle loop representing the bottom base of
            # the outlet.
            if isinstance(outlet, Cylinder):
                p_bottom_outlet, c_bottom_outlet, _ = gmsh_utils.add_circle_arc(
                    x=outlet.top_base_position[0],
                    y=outlet.top_base_position[1],
                    z=-outlet.height,
                    r=outlet.radius,
                )
            elif isinstance(outlet, Cube):
                p_bottom_outlet, c_bottom_outlet, _ = \
                    gmsh_utils.add_z_rectangle_loop(
                        x=outlet.top_base_position[0],
                        y=outlet.top_base_position[1],
                        z=-outlet.dimensions[2],
                        lx=outlet.dimensions[0],
                        ly=outlet.dimensions[1],
                    )

            # Add the walls of the outlet (cylindrical/cubic) block.
            gmsh_utils.add_cylinder_walls(p_bottom_outlet, c_bottom_outlet,
                                          p_top_outlet, c_top_outlet)

            # Add the loop representing the top base of the outlet to the list
            # of loops representing holes in the bottom base of the tank
            # cylinder.
            tank_base_hole_loops.append(l_top_outlet)

        # Add the top and bottom bases of the tank block, setting the loop
        # representing the top base of the outlet as a hole.
        if isinstance(shape, Cylinder):
            p_top, c_top, _, _ = gmsh_utils.add_circle(
                x=0,
                y=0,
                z=shape.height,
                r=shape.radius,
                hole_loops=[],
            )
            p_bottom, c_bottom, _, _ = gmsh_utils.add_circle(
                x=0,
                y=0,
                z=0,
                r=shape.radius,
                hole_loops=tank_base_hole_loops,
            )

        elif isinstance(shape, Cube):
            p_top, c_top, _, _ = gmsh_utils.add_z_rectangle(
                x=-shape.dimensions[0] / 2,
                y=-shape.dimensions[1] / 2,
                z=shape.dimensions[2],
                lx=shape.dimensions[0],
                ly=shape.dimensions[1],
                hole_loops=[],
            )
            p_bottom, c_bottom, _, _ = gmsh_utils.add_z_rectangle(
                x=-shape.dimensions[0] / 2,
                y=-shape.dimensions[1] / 2,
                z=0,
                lx=shape.dimensions[0],
                ly=shape.dimensions[1],
                hole_loops=tank_base_hole_loops,
            )

        # Add the walls of the tank (cylindrical/cubic) block.
        gmsh_utils.add_cylinder_walls(p_bottom, c_bottom, p_top, c_top)

    # Convert the msh file generated by gmsh to obj format.
    gmsh_utils.convert_msh_to_obj_file(path)


def _create_fluid_mesh_file(shape, fluid_level, margin, path: str):
    """Creates a mesh file for the fluid.
    
    The fluid is represented by a block with the same shape as the tank, but
    with a smaller height.

    Args:
        shape: Shape of the tank.
        fluid_level: Height of the fluid.
        margin: Margin to be added to the fluid block.
        path: Path of the file to be created.
    """

    if isinstance(shape, Cube):
        with gmsh_utils.gmshAPIWrapper():
            gmsh_utils.add_box(
                -shape.dimensions[0] / 2 + margin,
                -shape.dimensions[1] / 2 + margin,
                margin,
                shape.dimensions[0] - 2 * margin,
                shape.dimensions[1] - 2 * margin,
                fluid_level - 2 * margin,
            )

    elif isinstance(shape, Cylinder):
        with gmsh_utils.gmshAPIWrapper():
            gmsh_utils.add_cylinder(
                0,
                0,
                0 + margin,
                shape.radius - margin,
                fluid_level - margin,
            )
    else:
        raise ValueError(f"Invalid fluid shape `{shape}`.")

    # Convert the msh file generated by gmsh to obj format.
    gmsh_utils.convert_msh_to_obj_file(path)
