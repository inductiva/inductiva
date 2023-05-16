"""GROMACS module of the API"""

import pathlib
from typing import Optional

from inductiva import types
from inductiva.simulation import Simulator


class GROMACS(Simulator):
    """Class to invoke a GROMACS energy minimization on the API."""

    @property
    def api_method_name(self) -> str:
        return "md.gromacs.run_simulation"

    def run(
        self,
        input_directory: types.Path,
        method_name: str,
        output_directory: Optional[types.Path] = None,
        **gromacs_flags: Optional[types.Args],
    ) -> pathlib.Path:
        """Run a specefied gromacs method.

        Args:
            input_dir: Path to the directory containing the input files.
            method_name: Method to run.
            output_dir: Path to the directory where the output files will be
                stored.
            gromacs_flags: Flags to pass to the gromacs CLI.
        """
        return super().run(input_directory,
                           output_dir = output_directory,
                           method = method_name,
                           user_flags=gromacs_flags
                           )
