"""COAWST simulator module of the API."""

import os
import re
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

    def _regex_exists_in_file(self, file_path: str, pattern: str) -> bool:
        """
        Check if a given regular expression exists in a file.
        
        :param file_path: Path to the file.
        :param pattern: Regular expression pattern to search for.
        :return: True if the pattern exists, False otherwise.
        """
        with open(file_path, 'r', encoding='utf-8') as file:
            for line in file:
                if re.search(pattern, line):
                    return True
        return False

    def _validate_build_script(self, build_script: str) -> None:
        """
        Validate the build script by checking if required environment variables are correctly set.

        :param build_script: Path to the build script file.
        :raises ValueError: If any required variable is missing or incorrect.
        """
        checks = {
            #Check if MY_ROOT_DIR points to the correct folder
            r"\bexport MY_ROOT_DIR\s*=\s*/workdir/output/artifacts/__COAWST\b":
                "Wrong value for MY_ROOT_DIR in your build_coawst_script.\n"
                "MY_ROOT_DIR needs to be /workdir/output/artifacts/__COAWST",
            #Check if USE_MPI is on
            r"\bexport USE_MPI\s*=\s*on\b":
                "Wrong value for USE_MPI in your build_coawst_script.\n"
                "USE_MPI needs to be on",
            #Check if which_MPI is set to openmpi
            r"^\s*export\s+which_MPI\s*=\s*openmpi\b":
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

        if n_vcpus > on.n_vcpus.total:
            raise ValueError(
                "The number of virtual cpus asked surpasses the"
                " available virtual cpus for the selected resource.")

        build_script = os.path.join(input_dir, build_coawst_script)

        if input_dir is not None:
            config_filename = os.path.join(input_dir, sim_config_filename)
            if not os.path.isfile(config_filename):
                raise ValueError("The provided sim_config_filename is not "
                                 "present in your input directory.")
            if not os.path.isfile(build_script):
                raise ValueError("The provided build_coawst_script is not "
                                 "present in your input directory.")

        self._validate_build_script(build_script=str(build_script))

        mpi_kwargs = {}
        mpi_kwargs["use_hwthread_cpus"] = use_hwthread
        if n_vcpus is not None:
            mpi_kwargs["np"] = n_vcpus
        mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        commands = [
            #Copy COAWST source code to our input dir
            "cp -r /opt/COAWST /workdir/output/artifacts/__COAWST",
            #Compile COAWST
            f"bash {build_coawst_script}",
            Command(f"{coawst_bin} {sim_config_filename}",
                    mpi_config=mpi_config),
            "rm -r  __COAWST",
        ]

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           remote_assets=remote_assets,
                           resubmit_on_preemption=resubmit_on_preemption,
                           **kwargs)
