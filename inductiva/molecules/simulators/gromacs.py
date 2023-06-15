"""GROMACS module of the API"""

import pathlib
from typing import Optional

from inductiva import types
from inductiva.simulation import Simulator, Command


class GROMACS(Simulator):
    """Class to invoke any GROMACS command on the API."""

    @property
    def api_method_name(self) -> str:
        return "md.gromacs.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        method_name: str,
        track_logs: bool = False,
        output_directory: Optional[types.Path] = None,
        **gromacs_flags: Optional[str],
    ) -> pathlib.Path:
        """Run a specified gromacs method.

        Args:
            input_dir: Path to the directory containing the input files.
            method_name: Method to run.
            output_dir: Path to the directory where the output files will be
                stored.
            gromacs_flags: Flags to pass to the gromacs CLI.
        """
        return super().run(input_dir,
                           output_dir=output_directory,
                           track_logs=track_logs,
                           method=method_name,
                           user_flags=gromacs_flags)

    def run_async(
        self,
        input_dir: types.Path,
        method_name: str,
        **gromacs_flags: Optional[str],
    ) -> str:
        """Run a specified gromacs method asynchronously.

        Args:
            input_dir: Path to the directory containing the input files.
            method_name: Method to run.
            gromacs_flags: Flags to pass to the gromacs CLI.
        """
        return super().run_async(input_dir,
                                 method=method_name,
                                 user_flags=gromacs_flags)


class GROMACSCommand(Command):
    """Class for construction of GROMACS commands."""

    def __init__(self, method_name: str, **gromacs_flags: Optional[str]):
        """Command constructor.

        Args:
            method_name: The name of the method to be executed in GROMAS.
            **gromacs_flags: Additional keyword arguments to be passed 
            to the simulation API method.
        """
        super().__init__(method_name=method_name, **gromacs_flags)



