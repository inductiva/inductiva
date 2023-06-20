"""GROMACS module of the API"""

import pathlib
from typing import Optional

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
        track_logs: bool = False,
        method_name: Optional[str] = None,
        pipeline: Optional[list] = None,
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
                           pipeline=pipeline,
                           method=method_name,
                           user_flags=gromacs_flags)

    def run_async(
        self,
        input_dir: types.Path,
        pipeline: Optional[list] = None,
        method_name: Optional[str] = None,
        **gromacs_flags: Optional[str],
    ) -> str:
        """Run a specified gromacs method asynchronously.

        Args:
            input_dir: Path to the directory containing the input files.
            method_name: Method to run.
            gromacs_flags: Flags to pass to the gromacs CLI.
        """
        return super().run_async(input_dir,
                                 pipeline=pipeline,
                                 method=method_name,
                                 user_flags=gromacs_flags)
