"""Utils to create the different types of boundary conditions."""

from typing import List

import json

from inductiva.structures import bcs


class BoundaryConditionsCase:
    """Boundary conditions for a given case.

    The boundary conditions for a given case are characterized by a list of 
      Dirichlet and Neumann boundary conditions.

    Attributes:
        bcs_list (list): The boundary conditions objects.
    """

    def __init__(self, bcs_list: List[bcs.BoundaryCondition]) -> None:
        """Initializes a BoundaryConditionsCase object."""
        self.bcs_list = bcs_list

    def write_to_json(self, json_path: str) -> None:
        """Write the boundary conditions to JSON file.

        Args:
            json_path (str): The JSON file path.
        """

        # Divide the bcs_obj_list into two lists based on their types
        dirichlet_bcs = [
            bc for bc in self.bcs_list if isinstance(bc, bcs.DirichletBC)
        ]
        neumann_bcs = [
            bc for bc in self.bcs_list if isinstance(bc, bcs.NeumannBC)
        ]

        # Dirichlet dictionary
        if dirichlet_bcs is not None:
            list_of_dirichlet_bc_dicts = []
            for dirichlet_bc in dirichlet_bcs:
                dirichlet_bc_dict = dirichlet_bc.to_dict()
                list_of_dirichlet_bc_dicts.append(dirichlet_bc_dict)
            dirichlet_dict = {"dirichlet": list_of_dirichlet_bc_dicts}

            # Merge dictionaries: boundary conditions dictionary
            bcs_dict = {**dirichlet_dict}

        # Neumann dictionary
        if neumann_bcs is not None:
            list_of_neumann_bcs_dicts = []
            for neumann_bc in neumann_bcs:
                neumann_bc_dict = neumann_bc.to_dict()
                list_of_neumann_bcs_dicts.append(neumann_bc_dict)
            neumann_dict = {"neumann": list_of_neumann_bcs_dicts}

            # Merge dictionaries: boundary conditions dictionary
            bcs_dict = {**bcs_dict, **neumann_dict}

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(bcs_dict, write_file, indent=4)
