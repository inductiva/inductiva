"""SplisHSPlasH simulator module of the API."""
from typing import List, Optional

from inductiva import types, tasks, simulators


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
            Task object representing the simulation task.
        """

        self._input_files_exist(input_dir=input_dir,
                                remote_assets=remote_assets,
                                sim_config_filename=sim_config_filename)

        if vtk_to_obj and (vtk_dir is None or out_dir is None or particle_radius is None):
            raise ValueError("When using `vtk_to_obj=True` `vtk_dir` and "
                             "`out_dir` and `particle_radius` need to be "
                             "defined.")

        commands = [
            "cp /SPlisHSPlasH_CPU/bin/SPHSimulator .",
            f"./SPHSimulator {sim_config_filename} --no-gui --output-dir .",
            "rm SPHSimulator"
        ]

        if vtk_to_obj:
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
                "{}.obj"
            )
            
        return super().run(
            input_dir,
            commands=commands,
            storage_dir=storage_dir,
            on=on,
            resubmit_on_preemption=resubmit_on_preemption,
            remote_assets=remote_assets,
            project=project,
            **kwargs,
        )
