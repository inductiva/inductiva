"""Utils to create the material."""

import json


class IsotropicLinearElasticMaterial:
    """Linear elastic isotropic material.

    Attributes:
        young_modulus: A float representing the Young's modulus.
        poisson_ratio: A float representing the Poisson's ratio.
    """

    def __init__(self, young_modulus: float, poisson_ratio: float) -> None:
        """Initializes a IsotropicLinearElasticMaterial object."""
        self.young_modulus = young_modulus
        self.poisson_ratio = poisson_ratio

    def to_dict(self) -> dict:
        """Convert Material properties to a dictionary.

        Returns:
            material_dict: A dictionary containing the Material properties.
        """
        young_modulus_dict = {"young_modulus": self.young_modulus}
        poisson_ratio_dict = {"poisson_ratio": self.poisson_ratio}

        material_dict = {**young_modulus_dict, **poisson_ratio_dict}
        return material_dict

    def write_to_json(self, json_path: str) -> None:
        """Writes the material to JSON file.

        Args:
            json_path: A string representing the JSON file path.
        """
        material_dict = self.to_dict()

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(material_dict, write_file, indent=4)
