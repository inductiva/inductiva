"""Class to run commands on Quantum Espresso."""
from typing import List, Optional
import logging

from inductiva import types, tasks, simulators


@simulators.simulator.mpi_enabled
class QuantumEspresso(simulators.Simulator):
    """Class to run commands on Quantum Espresso."""

    COMMANDS = [
        "alpha2f", "dvscf_q2r", "head", "matdyn", "plan_avg", "pw", "rism1d",
        "turbo_spectrum", "average", "dynmat", "hp", "molecularnexafs",
        "plotband", "pw2bgw", "scan_ibrav", "upfconv", "band_interpolation",
        "epa", "ibrav2cell", "molecularpdos", "plotproj", "pw2critic", "simple",
        "virtual_v2", "bands", "epsilon", "initial_state", "neb", "plotrho",
        "pw2gt", "simple_bse", "wannier90", "bse_main", "ev", "kcw",
        "open_grid", "pmw", "pw2gw", "simple_ip", "wannier_ham", "casino2upf",
        "fermi_proj", "kcwpp_interp", "oscdft_et", "postahc", "pw2wannier90",
        "spectra_correction", "wannier_plot", "cell2ibrav", "fermi_velocity",
        "kcwpp_sh", "oscdft_pp", "postw90", "pw4gww", "sumpdos", "wfck2r", "cp",
        "fqha", "kpoints", "path_interpolation", "pp", "pwcond",
        "turbo_davidson", "wfdd", "cppp", "fs", "lambda", "pawplot", "ppacf",
        "pwi2xsf", "turbo_eels", "xspectra", "d3hess", "gww", "ld1", "ph",
        "pprism", "q2qstar", "turbo_lanczos", "dos", "gww_fit", "manycp",
        "phcg", "projwfc", "q2r", "turbo_magnon"
    ]

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the Quantum Espresso simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "quantumespresso"

    @property
    def name(self):
        """Get the name of the this simulator."""
        return "Quantum-Espresso"

    def available_commands(self) -> List[str]:
        """Get the list of available commands for Quantum Espresso.
        Returns:
            List of available commands for Quantum Espresso.
        """
        logging.info("\nAll available commands for Quantum Espresso have two"
                     " versions.\nMPI version and OpenMP version.\n"
                     "For the OpenMP version, append `_openmp` to the command "
                     "names (e.g., pw_openmp.x).\n"
                     "For the MPI version, just use the normal command names "
                     " (e.g., pw.x).")
        return self.COMMANDS

    def run(self,
            input_dir: Optional[str],
            commands: List[str],
            *,
            storage_dir: Optional[str] = "",
            on: types.ComputationalResources,
            extra_metadata: Optional[dict] = None,
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            project: Optional[str] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run.
            on: The computatißonal resource to launch the simulation on.
            storage_dir: Parent directory for storing simulation
                               results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
            project: Name of the project to which the task will be
                assigned. If None, the task will be assigned to
                the default project. If the project does not exist, it will be
                created.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata,
                           container_image=self._image_uri,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           **kwargs)
