"""COAWST simulator module of the API."""

import os
from typing import List, Optional

from inductiva import types, tasks, simulators
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class COAWST(simulators.Simulator):
    """Class to invoke a generic COAWST simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the COAWST simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "coawst"

    def _validate_build_script(self, build_script: str) -> None:
        """
        Validates the build script by checking if required environment variables
        are correctly set.

        :param build_script: Path to the build script file.
        :raises ValueError: If any required variable is missing or incorrect.
        """
        checks = {
            #Check if MY_ROOT_DIR points to the correct folder
            r"export\s+MY_ROOT_DIR\s*=\s*/workdir/output/artifacts/__COAWST":
                "Wrong value for MY_ROOT_DIR in your build_coawst_script.\n"
                "MY_ROOT_DIR needs to be /workdir/output/artifacts/__COAWST",
            #Check if USE_MPI is on
            r"export\s+USE_MPI\s*=\s*on":
                "Wrong value for USE_MPI in your build_coawst_script.\n"
                "USE_MPI needs to be on",
            #Check if which_MPI is set to openmpi
            r"export\s+which_MPI\s*=\s*openmpi":
                "Wrong value for which_MPI in your build_coawst_script.\n"
                "which_MPI needs to be openmpi"
        }

        for regex, error_msg in checks.items():
            if not self._regex_exists_in_file(build_script, regex):
                raise ValueError(error_msg)

    def run(self,
            input_dir: Optional[str],
            sim_config_filename: str,
            build_coawst_script: str,
            *,
            on: types.ComputationalResources,
            coawst_bin: str = "coawstM",
            init_commands: Optional[List[str]] = None,
            use_hwthread: bool = True,
            n_vcpus: Optional[int] = None,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
            sim_config_filename: Name of the simulation configuration file.
            build_coawst_script: Script used to build the COAWST executable.
            coawst_bin: Name of the COAWST binary to execute (coawstM by
            default).
            init_commands: List of helper commands to prepare things for your
                simulation. Used to copy files to and from
                `/workdir/output/artifacts/__COAWST`. It can also be used
                to run any helper function present within COAWST, like
                `scrip_coawst`.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
            (default), all vCPUs will be used.
            use_hwthread: If specified Open MPI will attempt to discover the
            number of hardware threads on the node, and use that as the
            number of slots available.
            other arguments: See the documentation of the base class.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
        """

        self._check_vcpus(n_vcpus, on)

        #only runs checks if we dont use remote assets
        if remote_assets is None:
            self._input_files_exist(input_dir=input_dir,
                                    sim_config_filename=sim_config_filename,
                                    build_coawst_script=build_coawst_script)

            build_script = os.path.join(input_dir, build_coawst_script)
            self._validate_build_script(build_script=str(build_script))

        mpi_kwargs = {}
        mpi_kwargs["use_hwthread_cpus"] = use_hwthread
        if n_vcpus is not None:
            mpi_kwargs["np"] = n_vcpus
        mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        # 34 selects dmpar for linux when compiling WRF
        compilation_command = Command(f"bash {build_coawst_script}", "34")

        commands = [
            #Copy COAWST source code to our input dir
            "cp -r /opt/COAWST /workdir/output/artifacts/__COAWST",
            "create_all_sim_links",
            #Compile COAWST
            compilation_command,
            Command(f"{coawst_bin} {sim_config_filename}",
                    mpi_config=mpi_config),
            "rm -r  __COAWST",
            "clean_all_sim_links"
        ]

        # Add init commands after building and before running the simulation
        if init_commands is not None:
            commands = commands[:1] + init_commands + commands[1:]

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           remote_assets=remote_assets,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
