"""GROMACS module of the API"""

from typing import List, Literal, Optional

from inductiva import types, tasks, simulators


class GROMACS(simulators.Simulator):
    """Class to invoke any GROMACS command on the API."""

    def __init__(self,
                 /,
                 version: Optional[str] = None,
                 use_dev: bool = False,
                 device: Literal["auto", "cpu", "gpu"] = "auto"):
        """Initialize the GROMACS simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
            device (str): Select between CPU or GPU for running the simulation.
                Default is "auto", which will auto-detect if the machine has a
                GPU and use if it available, otherwise use the CPU.
        """
        super().__init__(version=version, use_dev=use_dev, device=device)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "gromacs"

    def run(
        self,
        input_dir: Optional[str],
        commands: types.Commands,
        *,
        on: types.ComputationalResources,
        storage_dir: Optional[str] = "",
        resubmit_on_preemption: bool = False,
        remote_assets: Optional[List[str]] = None,
        project: Optional[str] = None,
        **kwargs,
    ) -> tasks.Task:
        """Run a list of GROMACS commands.

        Args:
            input_dir: Path to the directory containing the input files.
            on: The computational resource to launch the simulation on.
            commands: List of commands to run using the GROMACS simulator.
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
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           **kwargs)
