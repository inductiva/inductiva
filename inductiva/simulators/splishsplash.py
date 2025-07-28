"""SplisHSPlasH simulator module of the API."""
from typing import Optional, Tuple, Union

from inductiva import simulators, tasks, types


class SplishSplash(simulators.Simulator):
    """Class to invoke a generic SPlisHSPlasH simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the SPlisHSplasH simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "splishsplash"

    def run(
        self,
        input_dir: Optional[str],
        sim_config_filename: str,
        *,
        on: types.ComputationalResources,
        storage_dir: Optional[str] = "",
        resubmit_on_preemption: bool = False,
        remote_assets: Optional[Union[str, list[str]]] = None,
        project: Optional[str] = None,
        time_to_live: Optional[str] = None,
        vtk_to_obj: Optional[bool] = False,
        vtk_to_obj_vtk_dir: Optional[str] = None,
        vtk_to_obj_vtk_prefix: Optional[str] = "PartFluid_",
        vtk_to_obj_out_dir: Optional[str] = None,
        vtk_to_obj_particle_radius: Optional[float] = None,
        vtk_to_obj_smoothing_length: Optional[float] = 2.0,
        vtk_to_obj_cube_size: Optional[float] = 1.0,
        vtk_to_obj_surface_threshold: Optional[float] = 0.6,
        gen_gif: Optional[bool] = False,
        gen_gif_cam_pos: Tuple[float, float, float] = (4.0, 1.0, 4.0),
        gen_gif_cam_fp: Tuple[float, float, float] = (0.0, 0.0, 0.0),
        on_finish_cleanup: Optional[Union[str, list[str]]] = None,
        **kwargs,
    ) -> tasks.Task:
        """Run the SPlisHSPlasH simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
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
            vtk_to_obj: Whether to convert the output VTK files to OBJ meshes
                using marching cubes.
            vtk_to_obj_vtk_dir: Directory containing VTK files to be converted.
            vtk_to_obj_out_dir: Directory where the generated OBJ files will be
                stored.
            vtk_to_obj_vtk_prefix: Prefix of the VTK files that will be
                converted.
                Default: PartFluid_ 
            vtk_to_obj_particle_radius: The particle radius of the input data.
            vtk_to_obj_smoothing_length: The smoothing length radius used for
                the SPH kernel, the kernel compact support radius will be twice
                the smoothing length (in multiplies of the particle radius).
                Default: 2.0
            vtk_to_obj_cube_size: The cube edge length used for marching cubes
                in multiplies of the particle radius, corresponds to the cell
                size of the implicit background grid.
                Default: 1.0
            vtk_to_obj_surface_threshold: The iso-surface threshold for the
                density, i.e. the normalized value of the reconstructed density
                level that indicates the fluid surface (in multiplies of the
                rest density).
                Default: 0.6
            gen_gif: Whether to generate an animated GIF from the
                simulation output.
            gen_gif_cam_pos: The position of the camera when generating the GIF.  
                Default: (4.0, 1.0, 4.0).
            gen_gif_cam_fp: The point in space the camera looks at (focus
                point).  
                Default: (0.0, 0.0, 0.0).
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
        Returns:
            Task object representing the simulation task.
        """

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                sim_config_filename=sim_config_filename)

        if vtk_to_obj and (vtk_to_obj_vtk_dir is None or
                           vtk_to_obj_particle_radius is None):
            raise ValueError("When using `vtk_to_obj=True`, "
                             "`vtk_to_obj_vtk_dir` and "
                             "`vtk_to_obj_particle_radius` need to be defined.")

        commands = [
            "cp /SPlisHSPlasH_CPU/bin/SPHSimulator .",
            f"./SPHSimulator {sim_config_filename} --no-gui --output-dir .",
            "rm SPHSimulator"
        ]

        if vtk_to_obj:

            # If out_dir is not provided, will save in the same directory as
            # the vtk files
            if vtk_to_obj_out_dir is None:
                vtk_to_obj_out_dir = vtk_to_obj_vtk_dir

            commands.append(
                "splashsurf reconstruct "
                f"{vtk_to_obj_vtk_dir}/{vtk_to_obj_vtk_prefix}"
                "{}.vtk "
                f"-r={vtk_to_obj_particle_radius} "
                f"-l={vtk_to_obj_smoothing_length} "
                f"-c={vtk_to_obj_cube_size} "
                f"-t={vtk_to_obj_surface_threshold} "
                "--subdomain-grid=on --mesh-cleanup=on "
                "--mesh-smoothing-weights=on --mesh-smoothing-iters=25 "
                "--normals=on --normals-smoothing-iters=10 "
                f"-o {vtk_to_obj_out_dir}/{vtk_to_obj_vtk_prefix}_surface"
                "{}.obj")
        if gen_gif:
            commands.append("python3 /home/scripts/gen_gif.py ./vtk res.gif "
                            f"--cam_pos {gen_gif_cam_pos[0]} "
                            f"{gen_gif_cam_pos[1]} "
                            f"{gen_gif_cam_pos[2]} "
                            f"--cam_fp {gen_gif_cam_fp[0]} "
                            f"{gen_gif_cam_fp[1]} "
                            f"{gen_gif_cam_fp[2]}")

        return super().run(
            input_dir,
            commands=commands,
            storage_dir=storage_dir,
            on=on,
            resubmit_on_preemption=resubmit_on_preemption,
            remote_assets=remote_assets,
            project=project,
            time_to_live=time_to_live,
            on_finish_cleanup=on_finish_cleanup,
            **kwargs,
        )
