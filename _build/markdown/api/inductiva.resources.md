# resources

Base class for machine groups.

### *class* ElasticMachineGroup(machine_type: str, zone: str | None = None, provider: ProviderType | str = 'GCP', threads_per_core: int = 2, data_disk_gb: int = 10, max_idle_time: timedelta | int = 3, auto_terminate_ts: datetime | None = None, auto_terminate_minutes: int | None = None, spot: bool = True, mpi_config: MPIConfig = None, mpi_version: str = '4.1.6', np: int = None, use_hwthread_cpus: bool = True, auto_resize_disk_max_gb: int | None = None, min_machines: int = 1, max_machines: int = 2)

Bases: `BaseMachineGroup`

Create an ElasticMachineGroup object.

An ElasticMachineGroup is a set of identical machines that can
automatically scale based on CPU load. The group starts with a
minimum number of machines and adjusts its size as needed scaling
to the maximum number of machines, ensuring both optimal performance
and cost efficiency.
Note: The machine group becomes active after calling the ‘start’ method,
and billing commences once the machines are initiated.

* **Parameters:**
  * **machine_type** – The type of GC machine to launch. Ex: “e2-standard-4”.
    Check [https://cloud.google.com/compute/docs/machine-resource](https://cloud.google.com/compute/docs/machine-resource) for
  * **types.** (*more information about machine*)
  * **zone** – The zone where the machines will be launched.
  * **provider** – The cloud provider of the machine group.
  * **threads_per_core** – The number of threads per core (1 or 2).
  * **data_disk_gb** – The size of the disk for user data (in GB).
  * **max_idle_time** – Time without executing any task, after which the
    resource will be terminated.
  * **auto_terminate_ts** – Moment in which the resource will be
    automatically terminated.
  * **auto_resize_disk_max_gb** – The maximum size in GB that the hard disk
    of the cloud VM can reach. If set, the disk will be
    automatically resized, during the execution of a task, when the
    free space falls below a certain threshold. This mechanism helps
    prevent “out of space” errors, that can occur when a task
    generates a quantity of output files that exceeds the size of
    the local storage. Increasing disk size during task execution
    increases the cost of local storage associated with the VM,
    therefore the user must set an upper limit to the disk size, to
    prevent uncontrolled costs. Once that limit is reached, the disk
    is no longer automatically resized, and if the task continues to
    output files, it will fail.
  * **auto_terminate_minutes** – Duration, in minutes, the MPICluster will be
    kept alive. After auto_terminate_minutes minutes the machine
    will be terminated. This time will start counting after calling
    this method.
  * **min_machines** – The minimum number of available machines. This is
    a quantity of machines that will be started initially and the
    minimum available machines, even in cases of low CPU load.
  * **max_machines** – The maximum number of machines a machine group
    can scale up to.
  * **spot** – Whether to use spot machines.

#### QUOTAS_EXCEEDED_SLEEP_SECONDS *= 60*

#### \_\_init_\_(machine_type: str, zone: str | None = None, provider: ProviderType | str = 'GCP', threads_per_core: int = 2, data_disk_gb: int = 10, max_idle_time: timedelta | int = 3, auto_terminate_ts: datetime | None = None, auto_terminate_minutes: int | None = None, spot: bool = True, mpi_config: MPIConfig = None, mpi_version: str = '4.1.6', np: int = None, use_hwthread_cpus: bool = True, auto_resize_disk_max_gb: int | None = None, min_machines: int = 1, max_machines: int = 2) → None

#### active_machines_to_str() → str

Returns a string representation of the
number of machines currently running.

#### allow_auto_start *= True*

#### auto_resize_disk_max_gb *: int | None* *= None*

#### auto_terminate_minutes *: int | None* *= None*

#### auto_terminate_ts *: datetime | None* *= None*

#### *property* available_vcpus

Returns the maximum number of vCPUs that can be used on a task.

On a machine group with 2 machines, each with 4 vCPUs, this will return
4.
On an elastic machine group, this will also return 4.
On an MPI cluster this will return the total number of vcpus because
we can run on the total number of vcpus.

#### can_start_resource() → bool

Check if the resource can be started.

This method checks if the resource can be started by checking
the available quotas and resource usage.

* **Returns:**
  True if the resource can be started, False otherwise.
* **Return type:**
  bool

#### create_time *= None*

#### data_disk_gb *: int* *= 10*

#### estimate_cloud_cost(verbose: bool = True)

Estimates a cost per hour of min and max machines in US dollars.

these are the estimted costs of having minimum and the
maximum number of machines up in the cloud. The final cost will vary
depending on the total usage of the machines.

#### *classmethod* from_api_response(resp: dict)

Creates a MachineGroup object from an API response.

#### get_available_mpi_slots(use_hwthread: bool | None = None) → int

#### get_mpi_config()

Get the MPI configuration for the cluster.

#### gpu_count() → bool

Returns the number of GPUs available in the resource.

#### has_gpu() → bool

Check if the machine group has a GPU.

#### *property* id

#### *property* idle_time *: timedelta*

Resource idle time in seconds.

#### machine_type *: str*

#### max_idle_time *: timedelta | int* *= 3*

#### max_machines *: int* *= 2*

#### min_machines *: int* *= 1*

#### mpi_config *: MPIConfig* *= None*

#### mpi_version *: str* *= '4.1.6'*

#### *property* n_vcpus

Returns the number of vCPUs available in the resource.

Returns a tuple with the total number of vCPUs and the number of vCPUs
per machine. For a machine group with 2 machines, each with 4 vCPUs,
this will return (8, 4).
In a case of an Elastic machine group, the total number of vCPUs is
the maximum number of vCPUs that can be used at the same time.

#### *property* name

#### np *: int* *= None*

#### num_machines *= 0*

#### provider *: ProviderType | str* *= 'GCP'*

#### quota_usage *= {}*

#### quota_usage_table_str(resource_usage_header: str) → str

#### set_mpi_config(mpi_config: MPIConfig | None = None, mpi_version: str = '4.1.6', np: int | None = None, use_hwthread_cpus: bool = True)

Set the MPI configuration for the cluster.
:param mpi_config: An MPIConfig object containing the MPI configuration.

> Will take precedence over other arguments if provided.
* **Parameters:**
  * **mpi_version** – The version of MPI to be used on the machines.
  * **np** – The number of processes to use for MPI commands.
  * **use_hwthread_cpus** – Whether to use hyperthreading or not.

#### short_name() → str

#### spot *: bool* *= True*

#### start(wait_for_quotas: bool = False, verbose: bool = True)

Starts a machine group.

* **Parameters:**
  **wait_for_quotas** – If True, the method will wait for quotas to
  become available before starting the resource.

#### *property* started

#### terminate(verbose: bool = True)

Terminates a machine group.

#### threads_per_core *: int* *= 2*

#### *property* total_ram_gb

#### use_hwthread_cpus *: bool* *= True*

#### zone *: str | None* *= None*

### *class* MPICluster(machine_type: str, zone: str | None = None, provider: ProviderType | str = 'GCP', threads_per_core: int = 2, data_disk_gb: int = 10, max_idle_time: timedelta | int = 3, auto_terminate_ts: datetime | None = None, auto_terminate_minutes: int | None = None, spot: bool = True, mpi_config: MPIConfig = None, mpi_version: str = '4.1.6', np: int = None, use_hwthread_cpus: bool = True, num_machines: int = 2)

Bases: `BaseMachineGroup`

Create a MPICluster object.

A MPI cluster is a collection of homogenous machines all working together on
: a common task given the configurations that are launched in Google Cloud.
  Note: The cluster will be available only after calling ‘start’ method.
  The billing will start only after the machines are started.
  <br/>
  Args:
  : machine_type: The type of GC machine to launch. Ex: “e2-standard-4”.
    : Check [https://cloud.google.com/compute/docs/machine-resource](https://cloud.google.com/compute/docs/machine-resource) for
      information about machine types.
    <br/>
    zone: The zone where the machines will be launched.
    provider: The cloud provider of the machine group.
    threads_per_core: The number of threads per core (1 or 2).
    data_disk_gb: The size of the disk for user data (in GB).
    max_idle_time: Time without executing any task, after which the
    <br/>
    > resource will be terminated.
    <br/>
    auto_terminate_minutes: Duration, in minutes, the MPICluster will be
    : kept alive. After auto_terminate_minutes minutes the machine
      will be terminated. This time will start counting after calling
      this method.
    <br/>
    auto_terminate_ts: Moment in which the resource will be
    : automatically terminated.
    <br/>
    num_machines: The number of virtual machines to launch.

#### QUOTAS_EXCEEDED_SLEEP_SECONDS *= 60*

#### \_\_init_\_(machine_type: str, zone: str | None = None, provider: ProviderType | str = 'GCP', threads_per_core: int = 2, data_disk_gb: int = 10, max_idle_time: timedelta | int = 3, auto_terminate_ts: datetime | None = None, auto_terminate_minutes: int | None = None, spot: bool = True, mpi_config: MPIConfig = None, mpi_version: str = '4.1.6', np: int = None, use_hwthread_cpus: bool = True, num_machines: int = 2) → None

#### active_machines_to_str() → str

Return the number of machines currently running.

#### allow_auto_start *= True*

#### auto_resize_disk_max_gb *= None*

#### auto_terminate_minutes *: int | None* *= None*

#### auto_terminate_ts *: datetime | None* *= None*

#### *property* available_vcpus

Returns the number of vCPUs available to the resource.

For a mpi cluster with 2 machines, each with 4 vCPUs, this will
return 8.

#### can_start_resource() → bool

Check if the resource can be started.

This method checks if the resource can be started by checking
the available quotas and resource usage.

* **Returns:**
  True if the resource can be started, False otherwise.
* **Return type:**
  bool

#### create_time *= None*

#### data_disk_gb *: int* *= 10*

#### estimate_cloud_cost(verbose: bool = True)

Estimates a cost per hour of min and max machines in US dollars.

these are the estimted costs of having minimum and the
maximum number of machines up in the cloud. The final cost will vary
depending on the total usage of the machines.

#### *classmethod* from_api_response(resp: dict)

Creates a MachineGroup object from an API response.

#### get_available_mpi_slots(use_hwthread: bool | None = None) → int

#### get_mpi_config()

Get the MPI configuration for the cluster.

#### gpu_count() → bool

Returns the number of GPUs available in the resource.

#### has_gpu() → bool

Check if the machine group has a GPU.

#### *property* id

#### *property* idle_time *: timedelta*

Resource idle time in seconds.

#### machine_type *: str*

#### max_idle_time *: timedelta | int* *= 3*

#### mpi_config *: MPIConfig* *= None*

#### mpi_version *: str* *= '4.1.6'*

#### *property* n_vcpus

Returns the number of vCPUs available in the resource.

Returns a tuple with the total number of vCPUs and the number of vCPUs
per machine. For a machine group with 2 machines, each with 4 vCPUs,
this will return (8, 4).
In a case of an Elastic machine group, the total number of vCPUs is
the maximum number of vCPUs that can be used at the same time.

#### *property* name

#### np *: int* *= None*

#### num_machines *: int* *= 2*

#### provider *: ProviderType | str* *= 'GCP'*

#### quota_usage *= {}*

#### quota_usage_table_str(resource_usage_header: str) → str

#### set_mpi_config(mpi_config: MPIConfig | None = None, mpi_version: str = '4.1.6', np: int | None = None, use_hwthread_cpus: bool = True)

Set the MPI configuration for the cluster.
:param mpi_config: An MPIConfig object containing the MPI configuration.

> Will take precedence over other arguments if provided.
* **Parameters:**
  * **mpi_version** – The version of MPI to be used on the machines.
  * **np** – The number of processes to use for MPI commands.
  * **use_hwthread_cpus** – Whether to use hyperthreading or not.

#### short_name() → str

#### spot *: bool* *= True*

#### start(wait_for_quotas: bool = False, verbose: bool = True)

Starts a machine group.

* **Parameters:**
  **wait_for_quotas** – If True, the method will wait for quotas to
  become available before starting the resource.

#### *property* started

#### terminate(verbose: bool = True)

Terminates a machine group.

#### threads_per_core *: int* *= 2*

#### *property* total_ram_gb

#### use_hwthread_cpus *: bool* *= True*

#### zone *: str | None* *= None*

### *class* MachineGroup(machine_type: str, zone: str | None = None, provider: ProviderType | str = 'GCP', threads_per_core: int = 2, data_disk_gb: int = 10, max_idle_time: timedelta | int = 3, auto_terminate_ts: datetime | None = None, auto_terminate_minutes: int | None = None, spot: bool = True, mpi_config: MPIConfig = None, mpi_version: str = '4.1.6', np: int = None, use_hwthread_cpus: bool = True, byoc: bool = False, mg_name: str | None = None, auto_resize_disk_max_gb: int | None = None, num_machines: int = 1)

Bases: `BaseMachineGroup`

Create a MachineGroup object.

A machine group is a collection of homogenous machines with given the
configurations that are launched in Google Cloud.
Note: The machine group will be available only after calling ‘start’ method.
The billing will start only after the machines are started.

* **Parameters:**
  * **machine_type** – The type of GC machine to launch. Ex: “e2-standard-4”.
    Check [https://cloud.google.com/compute/docs/machine-resource](https://cloud.google.com/compute/docs/machine-resource) for
    information about machine types.
  * **zone** – The zone where the machines will be launched.
  * **provider** – The cloud provider of the machine group.
  * **threads_per_core** – The number of threads per core (1 or 2).
  * **data_disk_gb** – The size of the disk for user data (in GB).
  * **max_idle_time** – Time without executing any task, after which the
    resource will be terminated.
  * **auto_terminate_ts** – Moment in which the resource will be
    automatically terminated.
  * **auto_terminate_minutes** – Duration, in minutes, the MPICluster will be
    kept alive. After auto_terminate_minutes minutes the machine
    will be terminated. This time will start counting after calling
    this method.
  * **auto_resize_disk_max_gb** – The maximum size in GB that the hard disk
    of the cloud VM can reach. If set, the disk will be
    automatically resized, during the execution of a task, when the
    free space falls below a certain threshold. This mechanism helps
    prevent “out of space” errors, that can occur when a task
    generates a quantity of output files that exceeds the size of
    the local storage. Increasing disk size during task execution
    increases the cost of local storage associated with the VM,
    therefore the user must set an upper limit to the disk size, to
    prevent uncontrolled costs. Once that limit is reached, the disk
    is no longer automatically resized, and if the task continues to
    output files, it will fail.
  * **num_machines** – The number of virtual machines to launch.
  * **spot** – Whether to use spot machines.

#### QUOTAS_EXCEEDED_SLEEP_SECONDS *= 60*

#### \_\_init_\_(machine_type: str, zone: str | None = None, provider: ProviderType | str = 'GCP', threads_per_core: int = 2, data_disk_gb: int = 10, max_idle_time: timedelta | int = 3, auto_terminate_ts: datetime | None = None, auto_terminate_minutes: int | None = None, spot: bool = True, mpi_config: MPIConfig = None, mpi_version: str = '4.1.6', np: int = None, use_hwthread_cpus: bool = True, byoc: bool = False, mg_name: str | None = None, auto_resize_disk_max_gb: int | None = None, num_machines: int = 1) → None

#### active_machines_to_str() → str

Return the number of machines currently running.

#### allow_auto_start *= True*

#### auto_resize_disk_max_gb *: int | None* *= None*

#### auto_terminate_minutes *: int | None* *= None*

#### auto_terminate_ts *: datetime | None* *= None*

#### *property* available_vcpus

Returns the maximum number of vCPUs that can be used on a task.

On a machine group with 2 machines, each with 4 vCPUs, this will return
4.
On an elastic machine group, this will also return 4.
On an MPI cluster this will return the total number of vcpus because
we can run on the total number of vcpus.

#### byoc *: bool* *= False*

#### can_start_resource() → bool

Check if the resource can be started.

This method checks if the resource can be started by checking
the available quotas and resource usage.

* **Returns:**
  True if the resource can be started, False otherwise.
* **Return type:**
  bool

#### create_time *= None*

#### data_disk_gb *: int* *= 10*

#### estimate_cloud_cost(verbose: bool = True)

Estimates a cost per hour of min and max machines in US dollars.

these are the estimted costs of having minimum and the
maximum number of machines up in the cloud. The final cost will vary
depending on the total usage of the machines.

#### *classmethod* from_api_response(resp: VMGroupConfig)

Creates a MachineGroup object from an API response.

#### get_available_mpi_slots(use_hwthread: bool | None = None) → int

#### get_mpi_config()

Get the MPI configuration for the cluster.

#### gpu_count() → bool

Returns the number of GPUs available in the resource.

#### has_gpu() → bool

Check if the machine group has a GPU.

#### *property* id

#### *property* idle_time *: timedelta*

Resource idle time in seconds.

#### machine_type *: str*

#### max_idle_time *: timedelta | int* *= 3*

#### mg_name *: str | None* *= None*

#### mpi_config *: MPIConfig* *= None*

#### mpi_version *: str* *= '4.1.6'*

#### *property* n_vcpus

Returns the number of vCPUs available in the resource.

Returns a tuple with the total number of vCPUs and the number of vCPUs
per machine. For a machine group with 2 machines, each with 4 vCPUs,
this will return (8, 4).
In a case of an Elastic machine group, the total number of vCPUs is
the maximum number of vCPUs that can be used at the same time.

#### *property* name

#### np *: int* *= None*

#### num_machines *: int* *= 1*

#### provider *: ProviderType | str* *= 'GCP'*

#### quota_usage *= {}*

#### quota_usage_table_str(resource_usage_header: str) → str

#### set_mpi_config(mpi_config: MPIConfig | None = None, mpi_version: str = '4.1.6', np: int | None = None, use_hwthread_cpus: bool = True)

Set the MPI configuration for the cluster.
:param mpi_config: An MPIConfig object containing the MPI configuration.

> Will take precedence over other arguments if provided.
* **Parameters:**
  * **mpi_version** – The version of MPI to be used on the machines.
  * **np** – The number of processes to use for MPI commands.
  * **use_hwthread_cpus** – Whether to use hyperthreading or not.

#### short_name() → str

#### spot *: bool* *= True*

#### start(wait_for_quotas: bool = False, verbose: bool = True)

Starts a machine group.

#### *property* started

#### terminate(verbose: bool = True)

Terminates a machine group.

#### threads_per_core *: int* *= 2*

#### *property* total_ram_gb

#### use_hwthread_cpus *: bool* *= True*

#### zone *: str | None* *= None*

### get()

Returns a list of ‘Resource’ objects.

### get_by_name(machine_name: str)

Returns the machine group corresponding to machine_name.

<a id="module-inductiva.resources.utils"></a>

Functions to manage or retrieve user resources.

### estimate_machine_cost(machine_type: str, spot: bool = False, zone: str = None)

Estimate the cloud cost of one machine per hour in US dollars.

* **Parameters:**
  * **machine_type** – The type of GC machine to launch. Ex: “e2-standard-4”.
    Check [https://cloud.google.com/compute/docs/machine-resource](https://cloud.google.com/compute/docs/machine-resource) for
    more information about machine types.
  * **zone** – The zone where the machines will be launched.
  * **spot** – Whether to use spot machines.

### get_available_machine_types(provider: str | ProviderType = ProviderType.GCP, machine_families: List[str] | None = None, machine_configs: List[str] | None = None, vcpus_range: Tuple[int, int] | None = None, memory_range: Tuple[int, int] | None = None, price_range: Tuple[float, float] | None = None, spot: bool | None = None) → List[MachineTypeInfo]

Get all available machine types from a specified provider,
allowing for multiple filtering parameters.

* **Parameters:**
  * **provider** (*Union* *[**str* *,* *ProviderType* *]*) – The cloud provider for which to list available machine types.
    Defaults to ProviderType.GCP.
  * **machine_families** (*Optional* *[**List* *[**str* *]* *]*) – A list of machine families to filter the available machine types
    (e.g., “c2”, “c2d”, “n4”, etc).
  * **machine_configs** (*Optional* *[**List* *[**str* *]* *]*) – A list of specific machine configurations to filter the results
    (e.g., “highcpu”, “highmem”, “standard”).
  * **vcpus_range** (*Optional* *[**Tuple* *[**int* *,* *int* *]* *]*) – A tuple defining the range of virtual CPUs (vCPUs) to filter the
    machine types.
  * **memory_range** (*Optional* *[**Tuple* *[**int* *,* *int* *]* *]*) – A tuple defining the range of memory (in MB) to filter the machine
    types.
  * **price_range** (*Optional* *[**Tuple* *[**float* *,* *float* *]* *]*) – A tuple defining the price range (in the respective currency) to
    filter the machine types.
  * **spot** (*Optional* *[**bool* *]*) – If set to True, filters for spot instances; if False, filters for
    on-demand instances.
* **Returns:**
  A list of available machine types that match the specified criteria.
* **Return type:**
  List[MachineTypeInfo]

::docsbannersmall
::
