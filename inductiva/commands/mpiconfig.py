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

    def _replace_underscores_in_keys(self, d):
        """Replace underscores with hyphens in dictionary
        Since we are using kwargs we cant use flags like use-hwthread-cpus
        so we need to use use_hwthread_cpus as the parameter and then replace
        the underscores with hyphens before sending the request to the API.
        """
        new_dict = {}
        for key, value in d.items():
            new_key = key.replace('_', '-')
            if isinstance(value, dict):
                new_dict[new_key] = self._replace_underscores_in_keys(value)
            else:
                new_dict[new_key] = value
        return new_dict

    def to_dict(self):
        """Convert this MPIConfig to a dictionary representation.
        This method will replace the "_" with "-".
        """
        return self._replace_underscores_in_keys(self.__dict__)
