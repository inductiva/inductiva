"""Physical properties of different fluid types."""


class WATER:
    """Sets the physical properties of water."""

    def __init__(self) -> None:
        """Initializes a WATER object."""

        self.density = 1e3
        self.kinematic_viscosity = 1e-6
