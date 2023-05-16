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
        track_logs: bool = False,
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
<<<<<<< HEAD
        return super().run(input_directory,
                           output_dir = output_directory,
                           track_logs = track_logs,
                           method = method_name,
                           user_flags=gromacs_flags
                           )
=======
        return super().run(input_dir,
                           output_dir=output_dir,
                           track_logs=track_logs,
                           input_filename=sim_config_filename,
                           protein_filename=protein_filename,
                           topology_filename=topology_filename)

    def run_async(
        self,
        input_dir: types.Path,
        sim_config_filename: str,
        protein_filename: str,
        topology_filename: str,
    ) -> str:
        """Run the simulation asynchronously.
        
        Args:
            sim_config_filename: Name of the .mdp file containing the
                energy minimization parameters.
            protein_filename: Name of the file containing the protein
                .gro file.
            topology_filename: Name of the file containing the topology
                .top file.
            """

        return super().run_async(input_dir,
                                 input_filename=sim_config_filename,
                                 protein_filename=protein_filename,
                                 topology_filename=topology_filename)
>>>>>>> main
