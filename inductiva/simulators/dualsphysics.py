"""DualSPHysics simulator module of the API."""

from typing import List, Literal, Optional

from inductiva import types, tasks, simulators


class DualSPHysics(simulators.Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    def __init__(self,
                 /,
                 version: Optional[str] = None,
                 use_dev: bool = False,
                 device: Literal["auto", "cpu", "gpu"] = "auto"):
        """Initialize the DualSPHysics simulator.

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
        self.simulator_name_alias = "dualsphysics"

    def run(
        self,
        input_dir: Optional[str],
        shell_script: str,
        *,
        on: types.ComputationalResources,
        storage_dir: Optional[str] = "",
        resubmit_on_preemption: bool = False,
        remote_assets: Optional[List[str]] = None,
        project: Optional[str] = None,
        vtk_to_obj: Optional[bool] = False,
        vtk_dir: Optional[str] = None,
        prefix: Optional[str] = "PartFluid_",
        out_dir: Optional[str] = None,
        particle_radius: Optional[float] = None,
        smoothing_length: Optional[float] = 2.0,
        cube_size: Optional[float] = 1.0,
        surface_threshold: Optional[float] = 0.6,
        **kwargs,
    ) -> tasks.Task:
        """Executes a DualSPHysics simulation.

        Args:
            input_dir: Directory with simulation input files.
            shell_script: Path to the shell script to run the simulation.
            on: The computational resource to launch the simulation on.
            storage_dir: Directory for storing results.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            project: Name of the project to which the task will be
                assigned. If None, the task will be assigned to
                the default project. If the project does not exist, it will be
                created.
            vtk_to_obj: Whether to convert the output VTK files to OBJ meshes
                using marching cubes.
            vtk_dir: Directory containing VTK files to be converted.
            out_dir: Directory where the generated OBJ files will be stored.
            prefix: Prefix of the VTK files that will be converted.
                Default: PartFluid_ 
            particle_radius: The particle radius of the input data.
            smoothing_length: The smoothing length radius used for the SPH
                kernel, the kernel compact support radius will be twice the
                smoothing length (in multiplies of the particle radius).
                Default: 2.0
            cube_size: The cube edge length used for marching cubes in
                multiplies of the particle radius, corresponds to the cell size
                of the implicit background grid.
                Default: 1.0
            surface_threshold: The iso-surface threshold for the density, i.e.
                the normalized value of the reconstructed density level that
                indicates the fluid surface (in multiplies of the rest density).
                Default: 0.6
        Returns:
            tasks.Task: An object representing the simulation task.
        """

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                shell_script=shell_script)

        if vtk_to_obj and (vtk_dir is None or
                           particle_radius is None):
            raise ValueError("When using `vtk_to_obj=True` `vtk_dir`, "
                             "`out_dir` and `particle_radius` need to be "
                             "defined.")

        commands = [f"bash {shell_script}"]

        if vtk_to_obj:
            
            # If out_dir is not provided, will save in the same directory as
            # the vtk files
            if out_dir is None:
                out_dir = vtk_dir

            commands.append(
                f"splashsurf reconstruct {vtk_dir}/{prefix}"
                "{}.vtk "
                f"-r={particle_radius} "
                f"-l={smoothing_length} "
                f"-c={cube_size} "
                f"-t={surface_threshold} "
                "--subdomain-grid=on --mesh-cleanup=on "
                "--mesh-smoothing-weights=on --mesh-smoothing-iters=25 "
                "--normals=on --normals-smoothing-iters=10 "
                f"-o {out_dir}/{prefix}_surface"
                "{}.obj")
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           **kwargs)
