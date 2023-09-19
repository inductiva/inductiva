# pylint: disable=unused-argument
"""Deformable plate scenario."""

import os
from typing import List, Optional

from functools import singledispatchmethod

from inductiva import tasks, resources
from inductiva.scenarios import Scenario
from inductiva.simulators import Simulator
from inductiva.simulators import FEniCSx
from inductiva.types import Path

from inductiva.structures import bcs, holes, plates, materials
from . import mesh_utils
from . import bcs_utils
from . import geometry_utils

GEOMETRY_FILENAME = "geometry.json"
MESH_FILENAME = "mesh.msh"
BCS_FILENAME = "bcs.json"
MATERIAL_FILENAME = "material.json"


class DeformablePlate(Scenario):
    """Plate linear elastic scenario.

    The plate linear elastic scenario is characterized by the plate, holes,
    boundary conditions, and material properties. The plate and holes together
    form the geometry, which is used to create the mesh for the simulation.
    """

    valid_simulators = [FEniCSx]

    def __init__(
        self,
        plate: plates.RectangularPlate,
        holes_list: List[holes.Hole],
        bcs_list: List[bcs.BoundaryCondition],
        material: materials.IsotropicLinearElasticMaterial,
    ):
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
                 simulator: Simulator = FEniCSx(),
                 machine_group: Optional[resources.MachineGroup] = None,
                 run_async: bool = False) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: The simulator to use for the simulation.
            machine_group: The machine group to use for the simulation.
            run_async: Whether to run the simulation asynchronously.
            mesh_filename: Mesh filename.
            bcs_filename: Boundary conditions filename.
            material_filename: Material filename.
        """
        simulator.override_api_method_prefix("deformable_plate")
        task = super().simulate(simulator,
                                machine_group=machine_group,
                                run_async=run_async,
                                mesh_filename=MESH_FILENAME,
                                bcs_filename=BCS_FILENAME,
                                material_filename=MATERIAL_FILENAME)
        return task

    @singledispatchmethod
    def create_input_files(self, simulator: Simulator):
        pass


@DeformablePlate.create_input_files.register
def _(self,
      simulator: FEniCSx,
      input_dir: Path) -> None:
    """Creates FEniCSx simulation input files."""

    # Geometry file
    geometry_path = os.path.join(input_dir, GEOMETRY_FILENAME)
    self.geometry.write_to_json(geometry_path)

    # Mesh file
    mesh = mesh_utils.GmshMesh(self.geometry)
    mesh_path = os.path.join(input_dir, MESH_FILENAME)
    mesh.write_to_msh(mesh_path)

    # BCs file
    bcs_path = os.path.join(input_dir, BCS_FILENAME)
    self.bcs_case.write_to_json(bcs_path)

    # Material file
    material_path = os.path.join(input_dir, MATERIAL_FILENAME)
    self.material.write_to_json(material_path)
