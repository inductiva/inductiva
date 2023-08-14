"""Utils to create the material."""

import json


class IsotropicLinearElasticMaterial:
    """Linear elastic isotropic material.

    Attributes:
        young_modulus (float): Young's modulus of the material.
        poisson_ratio (float): Poisson's ratio of the material.
    """

    def __init__(self, young_modulus: float, poisson_ratio: float) -> None:
        """Initializes a IsotropicLinearElasticMaterial object."""
        self.young_modulus = young_modulus
        self.poisson_ratio = poisson_ratio

    def to_dict(self) -> dict:
        """Convert material properties to a dictionary.

        Returns:
            dict: Material properties.
        """
        return {
            "young_modulus": self.young_modulus,
            "poisson_ratio": self.poisson_ratio
        }

    def write_to_json(self, json_path: str) -> None:
        """Writes the material to JSON file.

        Args:
            json_path (str): JSON file path.
        """
        material_dict = self.to_dict()

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(material_dict, write_file, indent=4)


# Define standard materials
STRUCTURAL_STEEL = IsotropicLinearElasticMaterial(young_modulus=200000,
                                                  poisson_ratio=0.3)
ALUMINUM_ALLOY = IsotropicLinearElasticMaterial(young_modulus=71000,
                                                poisson_ratio=0.33)
GRAY_CAST_IRON = IsotropicLinearElasticMaterial(young_modulus=110000,
                                                poisson_ratio=0.28)
MAGNESIUM_ALLOY = IsotropicLinearElasticMaterial(young_modulus=45000,
                                                 poisson_ratio=0.35)
TITANIUM_ALLOY = IsotropicLinearElasticMaterial(young_modulus=96000,
                                                poisson_ratio=0.36)
COPPER_ALLOY = IsotropicLinearElasticMaterial(young_modulus=110000,
                                              poisson_ratio=0.34)
