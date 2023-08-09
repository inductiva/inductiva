"""Utils to create the different types of boundary conditions."""

from typing import Optional
from typing import List
from typing import Literal
from typing import Union

import json


class DirichletBC:
    """Dirichlet boundary condition.

    Attributes:
        boundary_name (str): The boundary name of the plate where the Dirichlet
          boundary condition will apply. 
          The available options are: left, top, right and bottom.
        displacement_x (float): The imposed displacement in x-direction. Only
          fill it if want to resprtit the displamente. 
          Fill this attribute if you want to specify the displacement.
        displacement_y (float): The imposed displacement in y-direction. Only
          fill it if want to resprtit the displamente.
          Fill this attribute if you want to specify the displacement.
    """

    def __init__(self,
                 boundary_name: Literal["left", "top", "right", "bottom"],
                 displacement_x: Optional[float] = None,
                 displacement_y: Optional[float] = None) -> None:
        """Initializes a DirichletBC object."""
        self.boundary_name = boundary_name
        self.displacement_x = displacement_x
        self.displacement_y = displacement_y


class NeumannBC:
    """Neumann boundary condition.

    Attributes:
        boundary_name (str): The boundary name of the plate where the Neumann 
          boundary condition will apply. 
          The available options are: left, top, right and bottom.
        tension_x (float): The tension value in x-direction.
        tension_y (float): The tension value in y-direction.
    """

    def __init__(self, boundary_name: Literal["left", "top", "right", "bottom"],
                 tension_x: float, tension_y: float) -> None:
        """Initializes a NeumannBC object."""
        self.boundary_name = boundary_name
        self.tension_x = tension_x
        self.tension_y = tension_y


class BoundaryConditionsCase:
    """Boundary conditions for a given case.

    The boundary conditions for a given case are characterized by a list of 
      Dirichlet and Tension boundary conditions.

    Attributes:
        bcs_objects_list (list): The boundary conditions objects.
    """

    def __init__(self, bcs_objects_list: List[Union[DirichletBC,
                                                    NeumannBC]]) -> None:
        """Initializes a BoundaryConditionsCase object."""
        self.bcs_objects_list = bcs_objects_list

    def write_to_json(self, json_path: str) -> None:
        """Write the boundary conditions to JSON file.

        Args:
            json_path (str): The JSON file path.
        """

        # Divide the bcs_obj_list into two lists based on their types
        dirichlet_objects_list = [
            obj for obj in self.bcs_objects_list
            if isinstance(obj, DirichletBC)
        ]
        neumann_objects_list = [
            obj for obj in self.bcs_objects_list if isinstance(obj, NeumannBC)
        ]

        # Dirichlet dictionary
        if dirichlet_objects_list is not None:
            dirichlet_dicts_list = []
            for dirichlet_object in dirichlet_objects_list:
                dirichlet_object_dict = vars(dirichlet_object)
                dirichlet_dicts_list.append(dirichlet_object_dict)
            dirichlet_dict = {"dirichlet": dirichlet_dicts_list}

            # Merge dictionaries: boundary conditions dictionary
            bcs_dict = {**dirichlet_dict}

        # Neumann dictionary
        if neumann_objects_list is not None:
            tension_dicts_list = []
            for neumann_object in neumann_objects_list:
                neumann_object_dict = vars(neumann_object)
                tension_dicts_list.append(neumann_object_dict)
            neumann_dict = {"neumann": tension_dicts_list}

            # Merge dictionaries: boundary conditions dictionary
            bcs_dict = {**bcs_dict, **neumann_dict}

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(bcs_dict, write_file, indent=4)
