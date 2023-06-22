"""GROMACS module of the API"""

import pathlib
from typing import Optional, List

from inductiva import types
from inductiva.simulation import Simulator


class GROMACS(Simulator):
    """Class to invoke any GROMACS command on the API."""

    @property
    def api_method_name(self) -> str:
        return "md.gromacs.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        commands: List[dict],
        track_logs: bool = False,
        output_dir: Optional[types.Path] = None,
    ) -> pathlib.Path:
        """Run a specified gromacs method.

        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run using the GROMACS simulatort.
            track_logs: If True, the logs of the remote execution will be
            logged. 
            output_dir: Path to the directory where the output files will be
            stored. If not provided, a timestamped directory will be created.
        """
        return super().run(input_dir,
                           output_dir=output_dir,
                           track_logs=track_logs,
                           commands=commands)

    def run_async(
        self,
        input_dir: types.Path,
        method_name: str = None,
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
