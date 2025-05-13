"""MOHID simulator module of the API."""
from typing import List, Optional

from inductiva import types, tasks, simulators

ALLOWD_MOHID_COMMANDS = ["MohidWater.exe", "MohidLand.exe"]


@simulators.simulator.mpi_enabled
class MOHID(simulators.Simulator):
    """Class to invoke a generic MOHID simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the MOHID simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "mohid"

    def run(self,
            input_dir: Optional[str],
            *,
            command: str = "MohidWater.exe",
            working_dir: str = "",
            on: types.ComputationalResources,
            run_ddc: bool = False,
            n_vcpus: Optional[int] = None,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            project: Optional[str] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
            command: MOHID command to run. Can be `MohidWater.exe` or
                `MohidLand.exe`.
            working_dir: Path (relative to the input directory) to the directory
                where the simulation will run. The command picked and
                domain consolidation will run on that folder.
            run_ddc: Boolean that will determine if we are going to run 
                domain consolidation (`MohidDDC.exe`) after the
                simulation ends. This will combine the results of the multiple
                MPI processes into one.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
            (default), all vCPUs will be used.
            other arguments: See the documentation of the base class.
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

        self._check_vcpus(n_vcpus, on)

        if command not in ALLOWD_MOHID_COMMANDS:
            raise ValueError(
                f"Command {command} not supported.\n"
                f"Available commands are: {ALLOWD_MOHID_COMMANDS}.")

        #We are using the intel MPI from oneapi
        #We don't have support for this MPI on the back end yet, this is why
        #we are not creating an MPICONFIG.
        mpi_command = "mpirun "
        if n_vcpus:
            mpi_command += f"-np {n_vcpus} "

        commands = [f"{mpi_command} /usr/local/bin/{command}"]

        if run_ddc:
            commands.append("/usr/local/bin/MohidDDC.exe")

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           remote_assets=remote_assets,
                           run_subprocess_dir=working_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           project=project,
                           **kwargs)
