# simulators

## AMR-Wind

AmrWind module of the API for numerical simulations of fluid flows.

### *class* AmrWind(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Bases: `Simulator`

Class to invoke a generic AmrWind simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Initialize the AmrWind simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.
  * **device** (*str*) – The device to use for the simulation. If “auto”,
    the device will be selected automatically. If “cpu”, the
    simulation will run on the CPU. If “gpu”, the simulation
    will run on the GPU.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the this simulator.

#### run(input_dir: str | None, , sim_config_filename: str, on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, use_hwthread: bool = True, n_vcpus: int | None = None, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **on** – The computational resource to launch the simulation on.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **storage_dir** – Path to the directory where the simulation results
    will be stored.
  * **resubmit_on_preemption** – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **\*\*kwargs** – Keyword arguments to be passed to the base class.

#### *property* version

Get the version of the simulator.

## CANS

CaNS module of the API for numerical simulations of fluid flows.

### *class* CaNS(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Bases: `Simulator`

Class to invoke a generic CaNS simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Initialize the CaNS simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.
  * **device** (*str*) – Select between CPU or GPU for running the simulation.
    Default is “auto”, which will auto-detect if the machine has a
    GPU and use if it available, otherwise use the CPU.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, sim_config_filename: str, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, use_hwthread: bool = True, n_vcpus: int | None = None, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## CM1

CM1 module of the API.

### *class* CM1(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic CM1 simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the CM1 simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, sim_config_filename: str, , base: str | None = None, init3d: str | None = None, init_terrain: str | None = None, init_surface: str | None = None, input_sounding: str | None = None, landuse: str | None = None, recompile: bool = False, on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, n_vcpus: int | None = None, use_hwthread: bool = True, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **base** – Path to the  base-state conditions file (base.F), if
    necessary. Will use the default base state file if not provided.
  * **init3d** – Path to the 3D initial conditions file (init3d.F), if you
    need to add perturbations to the base state. Will use the
    default initial conditions file if not provided.
  * **init_terrain** – Path to the terrain file. If you are
    using terrain, you will have to specify the terrain via the zs
    array in the file “init_terrain.F”. Will use the default terrain
    file if not provided.
  * **init_surface** – Path to the surface conditions file. If you are using
    surface fluxes of heat/moisture/momentum, then you might have
    to specify the horizontal distribution of several variables in
    the file init_surface.F. Will use the default surface
    conditions file if not provided.
  * **input_sounding** – Path to the sounding file. Used if you are supplying
    an external sounding file.
  * **landuse** – Path to the landuse file. If you are using surface fluxes
    of heat/momentum/moisture, or if you are using the atmospheric
    radiation scheme, then you need to specify the surface
    conditions. Used if you are supplying an external landuse file.
  * **recompile** – Recompile flag, disabled by default. The simulation will
    use the pre-compiled version of the CM1 code. The value is
    ignored if any of the files base, init3d, init_terrain,
    or init_surface is provided, as the simulator will always
    recompile the code in that case. If enabled, the simulator
    will recompile for the machine architecture of the
    computational resource, which may result in a performance
    improvement (usually not significant).
  * **on** – The computational resource to launch the simulation on.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified, Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **storage_dir** – Directory for storing simulation results.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## COAWST

COAWST simulator module of the API.

### *class* COAWST(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic COAWST simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the COAWST simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, sim_config_filename: str, , build_coawst_script: str | None = None, on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, coawst_bin: str = 'coawstM', init_commands: List[str] | None = None, cleanup_commands: List[str] | None = None, compile_simulator: bool = True, use_hwthread: bool = True, n_vcpus: int | None = None, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **build_coawst_script** – Script used to build the COAWST executable.
  * **coawst_bin** – str
    Name of the COAWST binary to execute (defaults to ‘coawstM’).
    This binary will be used whether the simulator is compiled or
    not. If compilation is skipped, the binary must
    already exist in the input_dir.
  * **init_commands** – List of helper commands to prepare things for your
    simulation. Used to copy files to and from
    /workdir/output/artifacts/_\_COAWST. It can also be used
    to run any helper function present within COAWST, like
    scrip_coawst.
    Will run before the compilation of the simulator.
  * **cleanup_commands** – List of helper commands to clean up things
    after your simulation finishes. Used to copy files from
    /workdir/output/artifacts/_\_COAWST to your input_dir. Or
    to delete unwanted files generated during the simulation.
    Will run after the simulation ends.
  * **compile_simulator** – If True, the simulator will be compiled using the provided
    build_coawst_script, and the simulation will run using the
    specified coawst_bin.
    If False, the simulation will be run directly using the
    precompiled coawst_bin binary, which should already exist in
    the input_dir with the same name as coawst_bin.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **arguments** (*other*) – See the documentation of the base class.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## CP2K

CP2K module of the API.

### *class* CP2K(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Bases: `Simulator`

Class to invoke a generic CP2K simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Initialize the CP2K simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.
  * **device** (*str*) – Select between CPU or GPU for running the simulation.
    Default is “auto”, which will auto-detect if the machine has a
    GPU and use if it available, otherwise use the CPU.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, sim_config_filename: str | None, n_vcpus: int | None = None, use_hwthread: bool = True, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **storage_dir** – Directory for storing simulation results.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## Delft3D

Delft3D module of the API.

### *class* Delft3D(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic Delft3D simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.
  * **device** (*str*) – Specifies whether to use the CPU or GPU version of
    the Docker image. If auto is picked we will pick the
    appropriate device based on the machine used to run the
    simulation, meaning, if the machine has a GPU we will pick
    GPU, otherwise we pick CPU.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , commands: List[str | inductiva.commands.Command] | None = None, shell_script: str | None = None, on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **storage_dir** – Directory for storing results.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **commands** – List of commands to run the simulation. Cannot be used
    with shell_script.
  * **shell_script** – Name of the shell script to run the simulation.
    Cannot be used with commands.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## DualSPHysics

DualSPHysics simulator module of the API.

### *class* DualSPHysics(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Bases: `Simulator`

Class to invoke a generic DualSPHysics simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Initialize the DualSPHysics simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.
  * **device** (*str*) – Select between CPU or GPU for running the simulation.
    Default is “auto”, which will auto-detect if the machine has a
    GPU and use if it available, otherwise use the CPU.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, shell_script: str, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, vtk_to_obj: bool | None = False, vtk_to_obj_vtk_dir: str | None = None, vtk_to_obj_vtk_prefix: str | None = 'PartFluid_', vtk_to_obj_out_dir: str | None = None, vtk_to_obj_particle_radius: float | None = None, vtk_to_obj_smoothing_length: float | None = 2.0, vtk_to_obj_cube_size: float | None = 1.0, vtk_to_obj_surface_threshold: float | None = 0.6, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Executes a DualSPHysics simulation.

* **Parameters:**
  * **input_dir** – Directory with simulation input files.
  * **shell_script** – Path to the shell script to run the simulation.
  * **on** – The computational resource to launch the simulation on.
  * **storage_dir** – Directory for storing results.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **vtk_to_obj** – Whether to convert the output VTK files to OBJ meshes
    using marching cubes.
  * **vtk_to_obj_vtk_dir** – Directory containing VTK files to be converted.
  * **vtk_to_obj_out_dir** – Directory where the generated OBJ files will be
    stored.
  * **vtk_to_obj_vtk_prefix** – Prefix of the VTK files that will be
    converted.
    Default: 

    ```
    PartFluid_
    ```
  * **vtk_to_obj_particle_radius** – The particle radius of the input data.
  * **vtk_to_obj_smoothing_length** – The smoothing length radius used for
    the SPH kernel, the kernel compact support radius will be twice
    the smoothing length (in multiplies of the particle radius).
    Default: 2.0
  * **vtk_to_obj_cube_size** – The cube edge length used for marching cubes
    in multiplies of the particle radius, corresponds to the cell
    size of the implicit background grid.
    Default: 1.0
  * **vtk_to_obj_surface_threshold** – The iso-surface threshold for the
    density, i.e. the normalized value of the reconstructed density
    level that indicates the fluid surface (in multiplies of the
    rest density).
    Default: 0.6
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
* **Returns:**
  An object representing the simulation task.
* **Return type:**
  tasks.Task

#### *property* version

Get the version of the simulator.

## FDS

FDS simulator module of the API.

### *class* FDS(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic FDS simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the FDS simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, sim_config_filename: str, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, smokeview_case_name: str | None = None, smokeview_video_sources: List[str] | None = None, n_vcpus: int | None = None, n_omp_threads: int | None = None, n_mpi_processes: int | None = None, use_hwthread: bool = True, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **smokeview_case_name** – 

    if provided, Smokeview will be automatically
    invoked after the simulation using the following command:
    > smokeview -runscript <smokeview_case_name>

    A .ssf (smokeview script) file with a matching name must be
    present in the input directory. The case name must also match
    the basename of the .smv file generated by the FDS
    simulation (without the .smv extension).
  * **smokeview_video_sources** – list of PNG filename prefixes to be used
    for generating videos from image sequences rendered by
    Smokeview via the RENDERALL instruction. For example, if
    Smokeview produces PNGs in the format “my_scenario_0001.png”,
    “my_scenario_0002.png”, etc., you should specify
    [“my_scenario”] to generate a video from that sequence.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **n_omp_threads** – Number of OpenMP threads to use in the simulation.
    If not provided, it defaults to the number of available vcpus
    of the machine group where the task will run.
  * **n_mpi_processes** – Number of MPI processes that will run the
    simulation. If not provided, it defaults to 1. Note that the
    number of MPI processes can’t exceed the number of meshes
    of the simulation case.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **arguments** (*other*) – See the documentation of the base class.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## FVCOM

FVCOM simulator module of the API.

### *class* FVCOM(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic FVCOM simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the FVCOM simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, debug: int = 0, model: str = '', case_name: str = '', use_hwthread: bool = True, create_namelist: str = '', n_vcpus: int | None = None, working_dir: str | None = '', storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **casename** – Name of the simulation case.
  * **on** – The computational resource to launch the simulation on. If None
    the simulation is submitted to a machine in the default pool.
  * **debug** – Debug level of the simulation (from 0 to 7).
  * **model** – At the current moment we provide users with two
  * **options** – 
    - None (default): Uses default fvcom binary.
    - ’estuary’: Uses the fvcom_estuary binary.

    The modules used to compile each binary can be found in the
    kutu repository, in the make.inc and make_estuary.inc files
    [https://github.com/inductiva/kutu/tree/main/simulators/fvcom/](https://github.com/inductiva/kutu/tree/main/simulators/fvcom/)
  * **create_namelist** – Used to create a namelist file for the simulation.
    Example: ‘create_namelist=hello’ will create hello_run.nml in
    the working_dir.
  * **working_dir** – Path (relative to the input directory) to the directory
    where the simulation nml file is located. If not provided, the
    input directory is used.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiates with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## GROMACS

GROMACS module of the API

### *class* GROMACS(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Bases: `Simulator`

Class to invoke any GROMACS command on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False, device: Literal['auto', 'cpu', 'gpu'] = 'auto')

Initialize the GROMACS simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.
  * **device** (*str*) – Select between CPU or GPU for running the simulation.
    Default is “auto”, which will auto-detect if the machine has a
    GPU and use if it available, otherwise use the CPU.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, commands: List[str | inductiva.commands.Command], , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run a list of GROMACS commands.

* **Parameters:**
  * **input_dir** – Path to the directory containing the input files.
  * **on** – The computational resource to launch the simulation on.
  * **commands** – List of commands to run using the GROMACS simulator.
  * **storage_dir** – Parent directory for storing simulation
    results.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## GX

GX module of the API for numerical simulations.

### *class* GX(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic GX simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the GX simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, sim_config_filename: str, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **sim_config_filename** – The name of the simulation configuration file.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## MOHID

MOHID simulator module of the API.

### *class* MOHID(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic MOHID simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the MOHID simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , command: str = 'MohidWater.exe', working_dir: str = '', on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, run_ddc: bool = False, n_vcpus: int | None = None, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **command** – MOHID command to run. Can be MohidWater.exe or
    MohidLand.exe.
  * **working_dir** – Path (relative to the input directory) to the directory
    where the simulation will run. The command picked and
    domain consolidation will run on that folder.
  * **run_ddc** – Boolean that will determine if we are going to run
    domain consolidation (MohidDDC.exe) after the
    simulation ends. This will combine the results of the multiple
    MPI processes into one.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## NWChem

NWChem simulator module of the API.

### *class* NWChem(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic NWChem simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the NWChem simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, sim_config_filename: str, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, n_vcpus: int | None = None, use_hwthread: bool = True, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## OpenFAST

OpenFAST module of the API for wind turbine simulations.

### *class* OpenFAST(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic OpenFAST simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the OpenFAST simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, commands: List[str | inductiva.commands.Command], , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **commands** – List of commands to run using the OpenFAST simulator.
  * **on** – The computational resource to launch the simulation on.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## OpenFOAM

OpenFOAM module of the API for fluid dynamics.

### *class* OpenFOAM(distribution: str = 'foundation', version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic OpenFOAM simulation on the API.

Users can choose between the ESI or the Foundation version
by selecting the version on the initiliasation. Be aware, that
some input files may only work for a specific version.

#### \_\_init_\_(distribution: str = 'foundation', version: str | None = None, use_dev: bool = False)

Initialize the OpenFOAM simulator.

* **Parameters:**
  * **distribution** (*str*) – The distribution of OpenFOAM to use. Available
    distributions are “foundation” and “esi”. Default is
    “foundation”.
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , commands: List[List[str | inductiva.commands.Command]] | None = None, shell_script: str | None = None, on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **commands** – List of commands to run using the OpenFOAM simulator.
  * **shell_script** – Path to a shell script (relative to input_dir) to run
    the simulation.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## OpenSees

OpenSees module of the API.

### *class* OpenSees(version: str | None = None, use_dev: bool = False, interface: Literal['python', 'tcl', 'eesd'] = 'python')

Bases: `Simulator`

Class to invoke a generic OpenSees simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False, interface: Literal['python', 'tcl', 'eesd'] = 'python')

Initialize the OpenSees simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.
  * **interface** (*str*) – The interface to use for interacting with the
    simulator. Can be either “python” (default) or “tcl”.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, sim_config_filename: str | None, n_vcpus: int | None = None, use_hwthread: bool = True, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **storage_dir** – Directory for storing simulation results.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## OpenTelemac

OpenTelemac module of the API.

### *class* OpenTelemac(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic OpenTelemac simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.
  * **device** (*str*) – Specifies whether to use the CPU or GPU version of
    the Docker image. If auto is picked we will pick the
    appropriate device based on the machine used to run the
    simulation, meaning, if the machine has a GPU we will pick
    GPU, otherwise we pick CPU.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, commands: List[str | inductiva.commands.Command], , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **storage_dir** – Directory for storing results.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **commands** – List of commands to run the simulation. Cannot be used
    with shell_script.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## Quantum Espresso

Class to run commands on Quantum Espresso.

### *class* QuantumEspresso(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to run commands on Quantum Espresso.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the Quantum Espresso simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### available_commands() → List[str]

Get the list of available commands for Quantum Espresso.
:returns: List of available commands for Quantum Espresso.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the this simulator.

#### run(input_dir: str | None, commands: List[str], , storage_dir: str | None = '', on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, extra_metadata: dict | None = None, resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.
:param input_dir: Path to the directory containing the input files.
:param commands: List of commands to run.
:param on: The computatißonal resource to launch the simulation on.
:param storage_dir: Parent directory for storing simulation

> results.
* **Parameters:**
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## Reef3D

Reef3D simulator module of the API.

### *class* REEF3D(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic REEF3D simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the REEF3D simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster | None, n_vcpus: int | None = None, use_hwthread: bool = True, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## SCHISM

SCHISM module of the API.

### *class* SCHISM(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic SCHISM simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the SCHISM simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, num_scribes: int = 1, storage_dir: str | None = '', use_hwthread: bool = True, n_vcpus: int | None = None, resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **num_scribes** – The num_scribes as per the simulator documentation.
  * **https** – //schism-dev.github.io/schism/master/getting-started/running-model.html
  * **storage_dir** – Directory for storing simulation results.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **n_vcpus** – Number of virtual cpus
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## SNL-SWAN

SNL SWAN module of the API.

### *class* SNLSWAN(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic SNL SWAN simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the SNL SWAN simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the this simulator.

#### run(input_dir: str | None, , sim_config_filename: str | None, remote_assets: str | list[str] | None = None, project: str | None = None, resubmit_on_preemption: bool = False, on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', n_vcpus: int | None = None, use_hwthread: bool = True, command: str = 'swanrun', time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **sim_config_filename** – Name of the simulation configuration file.
    Mandatory when using ‘swanrun’ command.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **on** – The computational resource to launch the simulation on.
  * **storage_dir** – Directory for storing simulation results.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **command** – The command to run the simulation. Default is ‘swanrun’.
    The user can also specify ‘swan.exe’.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## SPlisHSPlasH

SplisHSPlasH simulator module of the API.

### *class* SplishSplash(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic SPlisHSPlasH simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the SPlisHSplasH simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, sim_config_filename: str, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, vtk_to_obj: bool | None = False, vtk_to_obj_vtk_dir: str | None = None, vtk_to_obj_vtk_prefix: str | None = 'PartFluid_', vtk_to_obj_out_dir: str | None = None, vtk_to_obj_particle_radius: float | None = None, vtk_to_obj_smoothing_length: float | None = 2.0, vtk_to_obj_cube_size: float | None = 1.0, vtk_to_obj_surface_threshold: float | None = 0.6, gen_gif: bool | None = False, gen_gif_cam_pos: Tuple[float, float, float] = (4.0, 1.0, 4.0), gen_gif_cam_fp: Tuple[float, float, float] = (0.0, 0.0, 0.0), on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the SPlisHSPlasH simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **on** – The computational resource to launch the simulation on. If None
    the simulation is submitted to a machine in the default pool.
  * **storage_dir** – Directory for storing simulation results.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **vtk_to_obj** – Whether to convert the output VTK files to OBJ meshes
    using marching cubes.
  * **vtk_to_obj_vtk_dir** – Directory containing VTK files to be converted.
  * **vtk_to_obj_out_dir** – Directory where the generated OBJ files will be
    stored.
  * **vtk_to_obj_vtk_prefix** – Prefix of the VTK files that will be
    converted.
    Default: 

    ```
    PartFluid_
    ```
  * **vtk_to_obj_particle_radius** – The particle radius of the input data.
  * **vtk_to_obj_smoothing_length** – The smoothing length radius used for
    the SPH kernel, the kernel compact support radius will be twice
    the smoothing length (in multiplies of the particle radius).
    Default: 2.0
  * **vtk_to_obj_cube_size** – The cube edge length used for marching cubes
    in multiplies of the particle radius, corresponds to the cell
    size of the implicit background grid.
    Default: 1.0
  * **vtk_to_obj_surface_threshold** – The iso-surface threshold for the
    density, i.e. the normalized value of the reconstructed density
    level that indicates the fluid surface (in multiplies of the
    rest density).
    Default: 0.6
  * **gen_gif** – Whether to generate an animated GIF from the
    simulation output.
  * **gen_gif_cam_pos** – The position of the camera when generating the GIF.
    Default: (4.0, 1.0, 4.0).
  * **gen_gif_cam_fp** – The point in space the camera looks at (focus
    point).
    Default: (0.0, 0.0, 0.0).
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
* **Returns:**
  Task object representing the simulation task.

#### *property* version

Get the version of the simulator.

## SWAN

SWAN module of the API.

### *class* SWAN(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic SWAN simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the SWAN simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , sim_config_filename: str | None, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, resubmit_on_preemption: bool = False, on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', n_vcpus: int | None = None, use_hwthread: bool = True, command: str = 'swanrun', on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **sim_config_filename** – Name of the simulation configuration file.
    Mandatory when using ‘swanrun’ command.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **on** – The computational resource to launch the simulation on.
  * **storage_dir** – Directory for storing simulation results.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **command** – The command to run the simulation. Default is ‘swanrun’.
    The user can also specify ‘swan.exe’.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

## SWASH

SWASH module of the API.

### *class* SWASH(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic SWASH simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the SWASH simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, sim_config_filename: str, , command: str = 'swashrun', use_hwthread: bool = True, n_vcpus: int | None = None, storage_dir: str | None = '', on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, gfortran_unbuffered_all: bool = True, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **sim_config_filename** – Name of the simulation configuration file.
  * **on** – The computational resource to launch the simulation on.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **command** – The command to run the simulation. Default is ‘swashrun’.
    The user can also specify ‘swash.exe’.
  * **storage_dir** – Directory for storing simulation results.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **gfortran_unbuffered_all** – If True, enables immediate flushing of
    output, potentially degrading performance. However, disabling
    this might lead to error files not being written.

#### *property* version

Get the version of the simulator.

## XBeach

XBeach module of the API.

### *class* XBeach(version: str | None = None, use_dev: bool = False)

Bases: `Simulator`

Class to invoke a generic XBeach simulation on the API.

#### \_\_init_\_(version: str | None = None, use_dev: bool = False)

Initialize the XBeach simulator.

* **Parameters:**
  * **version** (*str*) – The version of the simulator to use. If None, the
    latest available version in the platform is used.
  * **use_dev** (*bool*) – Request use of the development version of
    the simulator. By default (False), the production version
    is used.

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, n_vcpus: int | None = None, use_hwthread: bool = True, sim_config_filename: str | None = None, export_vtk: bool = False, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.

* **Parameters:**
  * **input_dir** – Path to the directory of the simulation input files.
  * **on** – The computational resource to launch the simulation on.
  * **sim_config_filename** – Name of the simulation configuration file.
    Deprecated: This parameter is no longer used and will be removed
    in a future version. Please use params.txt in the input
    directory instead.
  * **n_vcpus** – Number of vCPUs to use in the simulation. If not provided
    (default), all vCPUs will be used.
  * **use_hwthread** – If specified Open MPI will attempt to discover the
    number of hardware threads on the node, and use that as the
    number of slots available.
  * **export_vtk** (*bool*) – If True, after the XBeach simulation finishes
    generate the VTK file via xbeach_animator –export-vtk.
  * **storage_dir** – Directory for storing simulation results.
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]
  * **arguments** (*other*) – See the documentation of the base class.

#### *property* version

Get the version of the simulator.

## Custom Image

It is also possible to use a custom image, from a public Docker registry.

<a id="module-inductiva.simulators.custom_image"></a>

Class to run commands on an custom image.

### *class* CustomImage(container_image: str)

Bases: `Simulator`

Class to run commands on an custom image.

#### \_\_init_\_(container_image: str)

Initialize the ArbitraryImage class.
Point to the API method to run a simulation.
:param container_image: The container image to use for the simulation.

> Example: container_image=”docker://inductiva/kutu:xbeach_v1.23”

#### *classmethod* get_supported_resources()

Get the supported computational resources for this simulator.

#### *property* image_uri

Get the image URI for this simulator.

#### *property* name

Get the name of the simulator.

#### run(input_dir: str | None, commands: List[str], , on: resources.MachineGroup | resources.ElasticMachineGroup | resources.MPICluster, storage_dir: str | None = '', resubmit_on_preemption: bool = False, remote_assets: str | list[str] | None = None, project: str | None = None, time_to_live: str | None = None, on_finish_cleanup: str | list[str] | None = None, \*\*kwargs) → [Task](inductiva.tasks.md#inductiva.tasks.task.Task)

Run the simulation.
:param input_dir: Path to the directory containing the input files.
:param commands: List of commands to run.
:param on: The computational resource to launch the simulation on.
:param storage_dir: Parent directory for storing simulation

> results.
* **Parameters:**
  * **resubmit_on_preemption** (*bool*) – Resubmit task for execution when
    previous execution attempts were preempted. Only applicable when
    using a preemptible resource, i.e., resource instantiated with
    spot=True.
  * **remote_assets** – Additional remote files that will be copied to
    the simulation directory.
  * **project** – Name of the project to which the task will be
    assigned. If None, the task will be assigned to
    the default project. If the project does not exist, it will be
    created.
  * **time_to_live** – Maximum allowed runtime for the task, specified as a
    string duration. Supports common time duration formats such as
    “10m”, “2 hours”, “1h30m”, or “90s”. The task will be
    automatically terminated if it exceeds this duration after
    starting.
  * **on_finish_cleanup** – 

    Optional cleanup script or list of shell commands to remove
    temporary or unwanted files generated during the simulation.
    This helps reduce storage usage by discarding unnecessary
    output.
    - If a string is provided, it is treated as the path to a shell
    script that must be included with the simulation files.
    - If a list of strings is provided, each item is treated as an
    individual shell command and will be executed sequentially.
    All cleanup actions are executed in the simulation’s working
    directory, after the simulation finishes.
    .. rubric:: Examples

    on_finish_cleanup = “my_cleanup.sh”

    on_finish_cleanup = [
    : “rm -rf temp_dir”,
      “rm -f logs/debug.log”

    ]

#### *property* version

Get the version of the simulator.

::docsbannersmall
::
