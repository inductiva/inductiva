"""WRF simulator module of the API."""

from typing import List, Optional, Union

from inductiva import simulators, tasks, types
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class WRF(simulators.Simulator):
    """Class to invoke a generic WRF simulation on the API."""

    VALID_CASE_NAMES = [
        "em_b_wave",
        "em_convrad",
        "em_esmf_exp",
        "em_fire",
        "em_grav2d_x",
        "em_heldsuarez",
        "em_hill2d_x",
        "em_les",
        "em_quarter_ss",
        "em_real",
        "em_scm_xy",
        "em_seabreeze2d_x",
        "em_squall2d_x",
        "em_squall2d_y",
        "em_tropical_cyclone",
    ]

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the WRF simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "wrf"

    def _validate_build_script(self, build_script: str) -> None:
        """
        Validates the build script by checking if required environment variables
        are correctly set.

        :param build_script: Path to the build script file.
        :raises ValueError: If any required variable is missing or incorrect.
        """

    def run(self,
            input_dir: Optional[str],
            *,
            on: types.ComputationalResources,
            use_hwthread: bool = True,
            case_name: str = "em_real",
            n_vcpus: Optional[int] = None,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[Union[str, list[str]]] = None,
            init_commands: Optional[List[str]] = None,
            gen_gif: bool = False,
            gen_gif_files: Optional[List[str]] = None,
            gen_gif_output_dir: Optional[str] = ".",
            gen_gif_fps: int = 3,
            gen_gif_variable: Optional[List[str]] = "RAINNC",
            project: Optional[str] = None,
            time_to_live: Optional[str] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
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
            init_commands: List of helper commands to prepare things for your
                simulation. Used to copy files from `/WRF`. It can also be used
                to run any helper function present within `/WRF/WPS`.
                Will run before the `wrf.exe` command.
            project: Name of the project to which the task will be
                assigned. If None, the task will be assigned to
                the default project. If the project does not exist, it will be
                created.
            gen_gif (bool): If True, generates a GIF from the simulation
                output files. Default is False.
            gen_gif_files (List[str]): List of ordered output files to be used
                for generating the GIF. Required if `gen_gif` is True.
                The order of the files will determine the order of frames in
                the GIF.
            gen_gif_output_dir (str): Directory where the GIF will be saved.
                Default is the current directory.
            gen_gif_fps (int): Frames per second for the GIF. Default is 3.
            gen_gif_variable (str): Variable to be used for generating the GIF.
                Default is "RAINNC".
            time_to_live: Maximum allowed runtime for the task, specified as a
                string duration. Supports common time duration formats such as
                "10m", "2 hours", "1h30m", or "90s". The task will be
                automatically terminated if it exceeds this duration after
                starting.
        """

        if case_name not in self.VALID_CASE_NAMES:
            raise ValueError(f"Invalid case name: {case_name}. "
                             f"Valid case names are: {self.VALID_CASE_NAMES}")

        if gen_gif and gen_gif_files is None:
            raise ValueError("If `gen_gif` is True, `gen_gif_files` must be "
                             "provided.")

        self._check_vcpus(n_vcpus, on)

        #only runs checks if we dont use remote assets

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                sim_config_filename="namelist.input")

        mpi_kwargs = {}
        mpi_kwargs["use_hwthread_cpus"] = use_hwthread
        if n_vcpus is not None:
            mpi_kwargs["np"] = n_vcpus
        mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        commands = [f"/scripts/create_links.sh /WRF/test/{case_name}"]

        if init_commands:
            commands += init_commands

        commands += [
            Command("./wrf.exe", mpi_config=mpi_config),
            "/scripts/delete_links.sh"
        ]

        if gen_gif:
            # if gen gif var is a single value turn it into a list
            if isinstance(gen_gif_variable, str):
                gen_gif_variable = [gen_gif_variable]
            for var in gen_gif_variable:
                files = " ".join(gen_gif_files)
                commands.append(
                    f"conda run -n wrf-env python /scripts/gen_gif.py "
                    f"--files {files} "
                    f"--output-dir {gen_gif_output_dir} "
                    f"--fps {gen_gif_fps} "
                    f"--var {var} ")

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           remote_assets=remote_assets,
                           resubmit_on_preemption=resubmit_on_preemption,
                           project=project,
                           time_to_live=time_to_live,
                           **kwargs)
