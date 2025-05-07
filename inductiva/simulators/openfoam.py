"""OpenFOAM module of the API for fluid dynamics."""
from typing import List, Optional

from inductiva import types, tasks, simulators

AVAILABLE_OPENFOAM_DISTRIBUTIONS = ["foundation", "esi"]


@simulators.simulator.mpi_enabled
class OpenFOAM(simulators.Simulator):
    """Class to invoke a generic OpenFOAM simulation on the API.

    Users can choose between the ESI or the Foundation version
    by selecting the version on the initiliasation. Be aware, that
    some input files may only work for a specific version.
    """

    def __init__(self,
                 /,
                 distribution: str = "foundation",
                 version: Optional[str] = None,
                 use_dev: bool = False):
        """Initialize the OpenFOAM simulator.

        Args:
            distribution (str): The distribution of OpenFOAM to use. Available
                distributions are "foundation" and "esi". Default is
                "foundation".
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        if distribution not in AVAILABLE_OPENFOAM_DISTRIBUTIONS:
            raise ValueError(
                f"Distribution '{distribution}' of OpenFOAM is not supported. "
                f"Available distributions are: "
                f"{AVAILABLE_OPENFOAM_DISTRIBUTIONS}")

        self._distribution = distribution

        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = f"openfoam_{distribution}"

    @property
    def name(self):
        """Get the name of the simulator."""
        return "OpenFOAM-" + self._distribution

    def run(self,
            input_dir: Optional[str],
            *,
            commands: Optional[List[types.Commands]] = None,
            shell_script: Optional[str] = None,
            on: types.ComputationalResources,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            project: Optional[str] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
            commands: List of commands to run using the OpenFOAM simulator.
            shell_script: Path to a shell script (relative to input_dir) to run
                the simulation.
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
            other arguments: See the documentation of the base class.
        """
        if not commands and not shell_script:
            raise ValueError("Either 'commands' or 'shell_script'"
                             " must be provided.")
        if commands and shell_script:
            raise ValueError("Only one of 'commands' or 'shell_script'"
                             " must be provided.")

        if shell_script:

            self._input_files_exist(input_dir=input_dir,
                                    remote_assets=remote_assets,
                                    shell_script=shell_script)
            commands = [f"bash {shell_script}"]

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           **kwargs)
