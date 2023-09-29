"""Utils to create the geometric case."""

from typing import List, Optional

import json

from inductiva.structures import holes, plates


class GeometricCase:
    """Geometric case.

    The geometric case is characterized by a plate and a set of holes.

    Attributes:
        plate (RectangularPlate): Rectangular plate object.
        holes_list (List[Hole]): The holes objects.
    """

    def __init__(self,
                 plate: plates.RectangularPlate,
                 holes_list: Optional[List[holes.Hole]] = None) -> None:
        """Initializes a GeometricCase object."""
        self.plate = plate
        self.holes_list = holes_list

    def write_to_json(self, json_path: str) -> None:
        """Writes the geometric case to JSON file.

        Args:
            json_path (str): The JSON file path.
        """

        # Gemeometric case dictionary
        plate_dict = {"plate": self.plate.to_dict()}
        holes_dict = {"holes": [hole.to_dict() for hole in self.holes_list]}
        geom_case_dictionary = {**plate_dict, **holes_dict}

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(geom_case_dictionary, write_file, indent=4)
