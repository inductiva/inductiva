# pylint: disable=unused-argument
"""Deformable plate scenario."""

import os
from typing import List, Optional

from inductiva import resources, scenarios, simulators, structures, tasks, types

from . import bcs_utils
from . import geometry_utils

GEOMETRY_FILENAME = "geometry.json"
BCS_FILENAME = "bcs.json"
MATERIAL_FILENAME = "material.json"


class DeformablePlate(scenarios.Scenario):
    """Plate linear elastic scenario.

    The plate linear elastic scenario is characterized by the plate, holes,
    boundary conditions, and material properties. The plate and holes together
    form the geometry, which is used to create the mesh for the simulation.
    """

    valid_simulators = [simulators.FEniCSx]

    def __init__(self, plate: structures.plates.RectangularPlate,
                 holes_list: List[structures.holes.Hole],
                 bcs_list: List[structures.bcs.BoundaryCondition],
                 material: structures.materials.IsotropicLinearElasticMaterial):
        """Initializes the plate linear elastic scenario.

        Args:
            plate (RectangularPlate): The rectangular plate geometry.
            holes_list (List[Hole]): Holes in the plate.
            bcs_list (List[BoundaryCondition]): The boundary conditions applied
              to the plate.
            material (IsotropicLinearElasticMaterial): The material properties
              of the plate.
            geometry (GeometricCase): The plate with holes geometry.
            bcs_case (BoundaryConditionsCase): The boudnary coinditions for the
             palte with holes.
        """
        self.plate = plate
        self.holes_list = holes_list
        self.bcs_list = bcs_list
        self.material = material
        self.geometry = geometry_utils.GeometricCase(plate=self.plate,
                                                     holes_list=self.holes_list)
        self.bcs_case = bcs_utils.BoundaryConditionsCase(bcs_list=self.bcs_list)

    def simulate(self,
                 simulator: simulators.Simulator = simulators.FEniCSx(),
                 machine_group: Optional[resources.MachineGroup] = None,
                 storage_dir: Optional[str] = "",
                 global_refinement_meshing_factor: float = 1.0,
                 local_refinement_meshing_factor: float = 0.0) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: The simulator to use for the simulation.
            machine_group: The machine group to use for the simulation.
            storage_dir: The parent directory where simulation
              results will be stored.
            global_refinement_meshing_factor (float): The refinement factor for
              global refinement of the mesh. A higher value results in a finer
              mesh overall, increasing the number of elements in the entire
              mesh, andvleading to a more detailed representation of the
              geometry. Use this factor when you want to globally refine the
              mesh uniformly, without specific local focus.
            local_refinement_meshing_factor (float): The refinement factor for
              local refinement of the mesh. This factor controls the local
              refinement level of the mesh and is typically used for refining
              specific regions or features of the mesh. A higher value for this
              factor indicates a finer mesh in the regions of interest,
              providing more detailed resolution around certain features. Use
              this factor when you want to focus on refining specific areas
              while keeping the rest of the mesh less refined.
        """
        simulator.override_api_method_prefix("deformable_plate")
        task = super().simulate(
            simulator,
            machine_group=machine_group,
            storage_dir=storage_dir,
            geometry_filename=GEOMETRY_FILENAME,
            bcs_filename=BCS_FILENAME,
            material_filename=MATERIAL_FILENAME,
            global_refinement_meshing_factor=global_refinement_meshing_factor,
            local_refinement_meshing_factor=local_refinement_meshing_factor)

        return task

    def create_input_files(self, simulator: simulators.FEniCSx,
                           input_dir: types.Path) -> None:
        """Creates FEniCSx simulation input files."""

        # Geometry file
        geometry_path = os.path.join(input_dir, GEOMETRY_FILENAME)
        self.geometry.write_to_json(geometry_path)

        # BCs file
        bcs_path = os.path.join(input_dir, BCS_FILENAME)
        self.bcs_case.write_to_json(bcs_path)

        # Material file
        material_path = os.path.join(input_dir, MATERIAL_FILENAME)
        self.material.write_to_json(material_path)
