"""FUNWAVE simulator module of the API."""

from typing import Optional, Union

from inductiva import simulators, tasks, types
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class FUNWAVE(simulators.Simulator):
    """Class to invoke a generic FUNWAVE simulation on the API."""

    # Saves in which line each flag is (in the Makefile)
    # Will be used to comment/uncomment the line if the flag is used or not
    FLAG_TO_LINE_MAP = {
        "COUPLING": 11,
        "ZALPHA": 12,
        "MANNING": 13,
        "VESSEL": 14,
        "METEO": 15,
        "WIND": 16,
        "SEDIMENT": 17,
        "CHECK_MASS_CONSERVATION": 18,
        "TMP": 19,
        "TRACKING": 20,
        "DEEP_DRAFT_VESSEL": 21,
    }

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the FUNWAVE simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "funwave"

    # pylint: disable=unused-argument
    def run(self,
            input_dir: Optional[str],
            sim_config_filename: str,
            *,
            COUPLING: Optional[bool] = False,
            ZALPHA: Optional[bool] = False,
            MANNING: Optional[bool] = False,
            VESSEL: Optional[bool] = False,
            METEO: Optional[bool] = False,
            WIND: Optional[bool] = False,
            SEDIMENT: Optional[bool] = False,
            CHECK_MASS_CONSERVATION: Optional[bool] = False,
            TRACKING: Optional[bool] = False,
            on: types.ComputationalResources,
            n_vcpus: Optional[int] = None,
            use_hwthread: bool = True,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[Union[str, list[str]]] = None,
            project: Optional[str] = None,
            time_to_live: Optional[str] = None,
            on_finish_cleanup: Optional[Union[str, list[str]]] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
            sim_config_filename: Name of the simulation configuration file.
            COUPLING : bool, optional
                Nesting mode.
            ZALPHA : bool, optional
                Activate Z-alpha parameterization for numerical dispersion and
                nonlinear wave handling.
            MANNING : bool, optional
                Use Manning formula for bottom friction.
            VESSEL : bool, optional
                Include shipwake module.
            METEO : bool, optional
                Include meteo tsunami module.
            WIND : bool, optional
                Include wind effect.
            SEDIMENT : bool, optional
                Include sediment and morphological module.
            CHECK_MASS_CONSERVATION : bool, optional
                Correct mass conservation problem caused by wetting/drying.
            TRACKING : bool, optional
                Include Lagrangian tracking module.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
                (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
                number of hardware threads on the node, and use that as the
                number of slots available.
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
            time_to_live: Maximum allowed runtime for the task, specified as a
                string duration. Supports common time duration formats such as
                "10m", "2 hours", "1h30m", or "90s". The task will be
                automatically terminated if it exceeds this duration after
                starting.
            on_finish_cleanup :
                Optional cleanup script or list of shell commands to remove
                temporary or unwanted files generated during the simulation.
                This helps reduce storage usage by discarding unnecessary
                output.
                - If a string is provided, it is treated as the path to a shell
                script that must be included with the simulation files.
                - If a list of strings is provided, each item is treated as an
                individual shell command and will be executed sequentially.
                All cleanup actions are executed in the simulation's working
                directory, after the simulation finishes.
                Examples:
                    on_finish_cleanup = "my_cleanup.sh"

                    on_finish_cleanup = [
                        "rm -rf temp_dir",
                        "rm -f logs/debug.log"
                    ]
        """

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                sim_config_filename=sim_config_filename)

        mpi_kwargs = {}
        mpi_kwargs["use_hwthread_cpus"] = use_hwthread
        if n_vcpus is not None:
            mpi_kwargs["np"] = n_vcpus

        mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        commands = ["cp /FUNWAVE-TVD-Version_3.6/Makefile ."]

        # Add sed commands for flags set to True
        for flag, line in self.FLAG_TO_LINE_MAP.items():
            if locals().get(flag, False):
                # If the flag is true it will remove the comment from the
                # Makefile according to the line in FLAG_TO_LINE_MAP
                commands.append(f"sed -i '{line}s/^# *//' Makefile")

        commands += [
            "make",
            Command(f"funwave-work/compiled_funwave {sim_config_filename}",
                    mpi_config=mpi_config),
            "rm -r funwave-work",
            "rm Makefile",
        ]

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           remote_assets=remote_assets,
                           resubmit_on_preemption=resubmit_on_preemption,
                           project=project,
                           time_to_live=time_to_live,
                           on_finish_cleanup=on_finish_cleanup,
                           **kwargs)
