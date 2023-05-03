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
        input_dir: types.Path,
        sim_config_filename: str,
        protein_filename: str,
        topology_filename: str,
        output_dir: Optional[types.Path] = None,
        track_logs: bool = False,
    ) -> pathlib.Path:
        """Run the energy minimization of a protein.

        Args:
            input_dir: Path to the directory containing the input files.
            sim_config_filename: Name of the .mdp file containing the
                energy minimization parameters.
            protein_filename: Name of the file containing the protein
                .gro file.
            topology_filename: Name of the file containing the topology
                .top file.
        """
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
