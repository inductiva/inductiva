"""CM1 module of the API."""
from typing import Optional, Union

from inductiva import simulators, tasks, types
from inductiva.commands.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


@simulators.simulator.mpi_enabled
class CM1(simulators.Simulator):
    """Class to invoke a generic CM1 simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the CM1 simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "cm1"

    def run(self,
            input_dir: Optional[str],
            sim_config_filename: str,
            *,
            base: Optional[str] = None,
            init3d: Optional[str] = None,
            init_terrain: Optional[str] = None,
            init_surface: Optional[str] = None,
            input_sounding: Optional[str] = None,
            landuse: Optional[str] = None,
            recompile: bool = False,
            on: types.ComputationalResources,
            n_vcpus: Optional[int] = None,
            use_hwthread: bool = True,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[Union[str, list[str]]] = None,
            project: Optional[str] = None,
            time_to_live: Optional[str] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            base: Path to the  base-state conditions file (base.F), if
                necessary. Will use the default base state file if not provided.
            init3d: Path to the 3D initial conditions file (init3d.F), if you
                need to add perturbations to the base state. Will use the
                default initial conditions file if not provided.
            init_terrain: Path to the terrain file. If you are
                using terrain, you will have to specify the terrain via the `zs`
                array in the file "init_terrain.F". Will use the default terrain
                file if not provided.
            init_surface: Path to the surface conditions file. If you are using
                surface fluxes of heat/moisture/momentum, then you might have
                to specify the horizontal distribution of several variables in
                the file `init_surface.F`. Will use the default surface
                conditions file if not provided.
            input_sounding: Path to the sounding file. Used if you are supplying
                an external sounding file.
            landuse: Path to the landuse file. If you are using surface fluxes 
                of heat/momentum/moisture, or if you are using the atmospheric
                radiation scheme, then you need to specify the surface
                conditions. Used if you are supplying an external landuse file.
            recompile: Recompile flag, disabled by default. The simulation will
                use the pre-compiled version of the CM1 code. The value is
                ignored if any of the files `base`, `init3d`, `init_terrain`,
                or `init_surface` is provided, as the simulator will always
                recompile the code in that case. If enabled, the simulator
                will recompile for the machine architecture of the
                computational resource, which may result in a performance
                improvement (usually not significant).
            on: The computational resource to launch the simulation on.
            sim_config_filename: Name of the simulation configuration file.
            n_vcpus: Number of vCPUs to use in the simulation. If not provided
                (default), all vCPUs will be used.
            use_hwthread: If specified, Open MPI will attempt to discover the
                number of hardware threads on the node, and use that as the
                number of slots available.
            storage_dir: Directory for storing simulation results.
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
            other arguments: See the documentation of the base class.
        """
        working_dir = "/workdir/output/artifacts"
        cm1_path = f"{working_dir}/__cm1"
        mpi_config = None

        self._check_vcpus(n_vcpus, on)

        files_to_check = {}
        # Check only files that are not None from:
        # base, init3d, init_terrain, init_surface, input_sounding, and
        # landuse
        for file in [
                base, init3d, init_terrain, init_surface, input_sounding,
                landuse, sim_config_filename
        ]:
            if file:
                files_to_check[file] = file
        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                **files_to_check)

        # Create MPI config
        mpi_kwargs = {"use_hwthread_cpus": use_hwthread}
        if n_vcpus is not None:
            mpi_kwargs["np"] = n_vcpus
        mpi_config = MPIConfig(version="4.1.6", **mpi_kwargs)

        # Setup commands for copying files and cleaning up
        copy_files_commands = []
        cleanup_commands = []

        if any((base, init3d, init_terrain, init_surface)):
            recompile = True

        if recompile:
            copy_files_commands.append(f"cp -r /cm1 {cm1_path}")
            cleanup_commands.append(f"rm -r {cm1_path}")

        # Copy the Fortan files updated by the user to the CM1 directory
        if base:
            copy_files_commands.append(f"cp -f {base} {cm1_path}/src")
        if init3d:
            copy_files_commands.append(f"cp -f {init3d} {cm1_path}/src")
        if init_terrain:
            copy_files_commands.append(f"cp -f {init_terrain} {cm1_path}/src")
        if init_surface:
            copy_files_commands.append(f"cp -f {init_surface} {cm1_path}/src")

        # Copy (and rename) the input_sounding and landuse files to the work dir
        # If they are already in the work dir with the correct name, we don't
        # need to copy
        if input_sounding and input_sounding != "input_sounding":
            copy_files_commands.append(
                f"cp -f {input_sounding} {working_dir}/input_sounding")
            cleanup_commands.append(f"rm {working_dir}/input_sounding")
        if landuse and landuse != "LANDUSE.TBL":
            copy_files_commands.append(
                f"cp -f {landuse} {working_dir}/LANDUSE.TBL")
            cleanup_commands.append(f"rm {working_dir}/LANDUSE.TBL")

        # Construct the commands to run the simulation
        bin_base_path = ""
        executable = "cm1.exe"
        commands = []
        commands.extend(copy_files_commands)

        if recompile:
            commands.append(f"make -C {cm1_path}/src")
            bin_base_path = f"{cm1_path}/run/"

        commands.append(
            Command(
                f"{bin_base_path}{executable} {sim_config_filename}",
                mpi_config=mpi_config,
            ))

        commands.extend(cleanup_commands)

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           time_to_live=time_to_live,
                           **kwargs)
