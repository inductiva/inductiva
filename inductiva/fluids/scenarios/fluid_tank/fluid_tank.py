"""Classes that define a fluid tank scenario and simulate it via API."""

from dataclasses import dataclass, field
import tempfile
from typing import List, Literal, Optional, Union
import shutil

from inductiva.types import Path
from inductiva.fluids.shapes import BaseShape
from inductiva.fluids.shapes import Rectangle
from inductiva.fluids.shapes import Circle
from inductiva.fluids.shapes import Cube
from inductiva.fluids.shapes import Cylinder
from inductiva.fluids.fluid_types import FluidType
from inductiva.fluids.fluid_types import WATER
from inductiva.fluids.simulators import SPlisHSPlasHParameters
from inductiva.fluids.simulators import DualSPHysicsParameters


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
        inlet: Optional[BaseTankInlet] = CircularTankInlet(radius=0.1,
                                                           position=[0, 0]),
        outlet: Optional[BaseTankOutlet] = CylindricalTankOutlet(
            radius=0.1, height=0.1, top_base_position=[0, 0]),
    ):
        self.shape = shape
        self.fluid = fluid
        self.fluid_level = fluid_level
        self.inlet = inlet
        self.outlet = outlet

    def simulate(
        self,
        device: Literal["cpu", "gpu"] = "cpu",
        engine: Literal["DualSPHysics", "SPlisHSPlasH"] = "SPlisHSPlasH",
        output_dir: Optional[Path] = None,
        engine_parameters: Union[
            DualSPHysicsParameters,
            SPlisHSPlasHParameters] = SPlisHSPlasHParameters,
    ):
        # Create a temporary directory to store simulation input files
        self.input_temp_dir = tempfile.TemporaryDirectory()  #pylint: disable=consider-using-with

        if engine.lower() == "splishsplash" and \
            isinstance(engine_parameters, SPlisHSPlasHParameters):
            sim_output_path = self._splishsplash_simulation()
        elif engine.lower() == "dualsphysics" and \
            isinstance(engine_parameters, DualSPHysicsParameters):
            sim_output_path = self._dualsphysics_simulation()
            raise NotImplementedError(
                f"Simulations of the fluid tank scenario with engine `{engine}` are not supported."
            )
        else:
            raise ValueError("Entered `engine` does not exist or it \
                             does not match with `engine_parameters` class")

        # Delete temporary input directory
        self.input_temp_dir.cleanup()

        return ""  # SimulationOutput(sim_output_path)

    def _splishsplash_simulation(self):
        """Runs SPlisHSPlasH simulation via API."""

        input_dir = self.input_temp_dir.name

        fluid_block_dir = os.path.dirname(__file__)
        unit_box_file_path = os.path.join(fluid_block_dir,
                                          UNIT_BOX_MESH_FILENAME)
        shutil.copy(unit_box_file_path, input_dir)

        fluid_margin = 2 * self.particle_radius

        replace_params_in_template_file(
            templates_dir=fluid_block_dir,
            template_filename=SPLISHSPLASH_TEMPLATE_FILENAME,
            params={
                "simulation_time": self.simulation_time,
                "time_step": TIME_STEP,
                "particle_radius": self.particle_radius,
                "data_export_rate": 1 / self.engine_parameters.output_time_step,
                "tank_filename": UNIT_BOX_MESH_FILENAME,
                "tank_dimensions": TANK_DIMENSIONS,
                "fluid_filename": UNIT_BOX_MESH_FILENAME,
                "fluid": self.fluid,
                "fluid_position": [fluid_margin] * 3,
                "fluid_dimensions": [
                    dimension - 2 * fluid_margin
                    for dimension in self.dimensions
                ],
                "fluid_velocity": self.initial_velocity,
            },
            output_file_path=os.path.join(input_dir,
                                          SPLISHSPLASH_INPUT_FILENAME),
        )

        logging.info("Estimated number of particles %d",
                     self.estimate_num_particles())
        logging.info("Estimated number of time steps %s",
                     math.ceil(self.simulation_time / TIME_STEP))
        logging.info(
            "Number of output time steps %s",
            math.ceil(self.simulation_time /
                      self.engine_parameters.output_time_step))

        sim_output_path = inductiva.sph.splishsplash.run_simulation(
            sim_dir=input_dir,
            input_filename=SPLISHSPLASH_INPUT_FILENAME,
            device=self.device,
            output_dir=self.output_dir)

        convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(sim_output_path, "vtk"),
            output_time_step=self.engine_parameters.output_time_step,
            netcdf_data_dir=os.path.join(sim_output_path, "netcdf"))

        return sim_output_path
