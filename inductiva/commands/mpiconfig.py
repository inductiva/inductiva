"""MPIConfig class"""


class MPIConfig:
    """MPIConfig class"""

    def __init__(self, version: str, **kwargs):
        """
        Args:
          version: The version of the MPI library to use.
            At the moment, only 1.10.7 and 4.1.6 are supported.
          **kwargs: Additional parameters or flags that Openmpi uses.
            parameters that have "-" should be typed with "_" instead.
            Example: np=2, use_hwthread_cpus=True
        """
        self.version = version
        self.options = kwargs

    def to_dict(self):
        """Convert this MPIConfig to a dictionary representation.
        This method will replace the "_" with "-".
        """
        return {
            "version": self.version,
            "options": {
                k.replace("_", "-"): v for k, v in self.options.items()
            }
        }

    @classmethod
    def from_dict(cls, d: dict) -> "MPIConfig":
        """
        Create an MPIConfig object from a dictionary.

        Args:
            d: A dictionary with keys "version" and "options".

        Returns:
            An MPIConfig object.

        Example:
        >>> d = {
        ...   "version": "4.1.6",
        ...   "options": {"np": 4, "use-hwthread-cpus": True}
        ... }
        >>> cfg = MPIConfig.from_dict(d)
        >>> cfg.version
        '4.1.6'
        >>> cfg.options
        {'np': 4, 'use_hwthread_cpus': True}
        """
        version = d.get("version")
        options = d.get("options", {})

        # Replace "-" back to "_"
        normalized_options = {k.replace("-", "_"): v for k, v in options.items()}

        return cls(version, **normalized_options)
