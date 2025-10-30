# tasks

## inductiva.tasks methods

Methods to interact with the tasks submitted to the API.

### get(last_n: int = 5, status: str | None = None, project: str | None = None) → List[[Task](#inductiva.tasks.task.Task)]

Get the last N tasks of a user.

This function fetches info about the last N tasks (with respect to
submission time) of a user to stdout, sorted by submission time with the
most recent first.
A status can be specified to filter to get only tasks with that status, in
which case the last N tasks with that status will be listed.
The number of tasks can be less than N if the aren’t enough tasks that match
the specified criteria.

Similar to the inductiva.task.list() function, but instead of printing
to stdout, returns a list of Task objects which can be used to perform
further operations on those tasks (e.g., download outputs,
kill unfinished tasks, etc.).

Example usage:
: # get the last 5 tasks that haven’t started yet and kill them
  tasks = inductiva.tasks.get(5, status=”submitted”)
  for task in tasks:
  <br/>
  > task.kill()

* **Parameters:**
  * **last_n** – The number of most recent tasks with respect to submission
    time to fetch. If filtering criteria (currently status is available)
    is specified, most recent N tasks that match that criteria will be
    listed. The actual number of tasks may be less if there
    aren’t enough tasks available.
  * **status** – The status of the tasks to get. If None, tasks with any status
    will be returned.
  * **project** – The project from which to fetch. If None, fetches from all
    projects.
* **Returns:**
  List of Task objects.

### get_all(status: str | None = None, project: str | None = None) → List[[Task](#inductiva.tasks.task.Task)]

Get all tasks of a user.

This function fetches all tasks of a user, sorted by submission
time with the most recent first. If status is specified, only
tasks with that status will be fetched.
:param status: The status of the tasks to get. If None, tasks with any status

> will be returned.
* **Returns:**
  List of dictionaries with information about the tasks.

### get_tasks(last_n: int = 10, project: str | None = None, status: str | None = None)

Get the last N submitted tasks.

Get the last N submitted tasks, eventually filtered by status.
By default, only the last 10 submitted tasks are returned,
irrespectively of their status.

* **Parameters:**
  * **last_n** (*int*) – The number of tasks with repect to the submission
    time to fectch. If last_n<=0 we fetch all tasks submitted
    to the project.
  * **status** – Status of the tasks to get. If None, tasks with any
    status will be returned.

### to_dict(list_of_tasks: Iterable[[Task](#inductiva.tasks.task.Task)]) → Mapping[str, List[Any]]

Converts an Iterable of tasks to a dictionary with all the
relevant information for all the tasks.

> Args:
> : list_of_tasks: An Iterable of tasks.

> Returns:
> : A dictionary with all the relevant information for
>   all the tasks. Example: { “ID”: [1, 2, 3],
>   “Simulator”: [“reef3d”, “reef3d”, “reef3d”], … }

## Tasks class

Manage running/completed tasks on the Inductiva API.

### *class* Task(task_id: str)

Bases: `object`

Represents a running/completed task on the Inductiva API.

Example usage:

```python
task = simulator.run(...)
task.wait()
info = task.get_info() # dictionary with info about the task
task.download_outputs(
    filenames=["file1.txt", "file2.dat"] # download only these files
)
```

#### \_\_init_\_(task_id: str)

Initialize the instance from a task ID.

#### *async* close_stream()

Close the stream to the task.

#### download_inputs(filenames: List[str] | None = None, input_dir: str | None = None, uncompress: bool = True, rm_downloaded_zip_archive: bool = True, rm_remote_files: bool = False) → Path | None

Download input files of the task.

* **Parameters:**
  * **filenames** – List of filenames to download. If None or empty, all
    files are downloaded.
  * **input_dir** – Directory where to download the files. If None, the
    files are downloaded to the default directory. The default is
    {inductiva.get_output_dir()}/{task_id}/inputs/.
  * **uncompress** – Whether to uncompress the archive after downloading it.
  * **rm_downloaded_zip_archive** – Whether to remove the archive after
    uncompressing it. If uncompress is False, this argument is
    ignored.
  * **rm_remote_files** – Whether to remove all task files from remote
    storage after the download is complete. Only used if filenames
    is None or empty (i.e., all input files are downloaded).

#### download_outputs(filenames: List[str] | None = None, output_dir: str | None = None, uncompress: bool = True, rm_downloaded_zip_archive: bool = True, rm_remote_files: bool = False) → Path | None

Download output files of the task.

* **Parameters:**
  * **filenames** – List of filenames to download. If None or empty, all
    files are downloaded.
  * **output_dir** – Directory where to download the files. If None, the
    files are downloaded to the default directory. The default is
    {inductiva.get_output_dir()}/{task_id}/outputs/.
  * **uncompress** – Whether to uncompress the archive after downloading it.
  * **rm_downloaded_zip_archive** – Whether to remove the archive after
    uncompressing it. If uncompress is False, this argument is
    ignored.
  * **rm_remote_files** – Whether to remove all task files from remote
    storage after the download is complete. Only used if filenames
    is None or empty (i.e., all output files are downloaded).

#### *classmethod* from_api_info(info: Task) → [Task](#inductiva.tasks.task.Task)

#### get_computation_time(cached: bool = False) → float | None

Get the time the computation of the task took to complete.

* **Returns:**
  The task computation time if the task is already started or in a
  terminal state, or None otherwise.

#### get_info() → TaskInfo

Get a dictionary with information about the task.

Information includes e.g., “task_id”, “status”, timestamps
(“create_time”, “input_submit_time, “start_time”, “end_time”),
among others.

This method issues a request to the API.

#### get_input_url() → str | None

Get a public URL to download the input files of the task.

* **Returns:**
  The URL to download the input files of the task, or None

#### get_machine_type() → str | None

Get the machine type used in the task.

Streamlines the process of obtaining the task info, extracting the
machine type from the comprehensive task info.

* **Returns:**
  The machine type, or None if a machine hasn’t been assigned yet.

#### get_metadata() → Dict[str, str]

Get the metadata associated with the task.

* **Returns:**
  A dictionary with the custom metadata previously set on this task.

#### get_output_info() → TaskOutputInfo

Get information about the output files of the task.

* **Returns:**
  An instance of the OutputInfo class, which can be used to
  access info about the output archive (number of files, total
  compressed size, total uncompressed size) and information about
  each file (name, size, compressed size). It can also be used to
  print that information in a formatted way.

#### get_output_url() → str | None

Get a public URL to download the output files of the task.

* **Returns:**
  The URL to download the output files of the task, or None
  if the

#### get_position_in_queue() → int | None

Get the position of the task in the queue.

This method issues a request to the API.

#### get_simulator_name() → str

#### get_status() → TaskStatusCode

Get status of the task.

This method issues a request to the API and updates the task info
is_terminal. The api call to get the status now returns the status and
is_terminated.

#### get_storage_path() → str

Get the path to this task’s directory in the user’s remote storage.

* **Returns:**
  String with the path to the task’s directory in remote storage.

#### get_total_time(cached: bool = False) → float | None

Get the total time the task workflow took to complete.

* **Returns:**
  The task total duration since it was created, or None if the
  metric is not available or can’t be computed.

#### *property* info *: TaskInfo*

Get information about the task.

It contains cached information about the task from the latest call to
get_info, therefore it can be outdated.

#### is_failed() → bool

Validate if the task has failed.

This method issues a request to the API.

#### is_running() → bool

Validate if the task is running.

This method issues a request to the API.

#### is_terminal() → bool

Check if the task is in a terminal status.

This method issues a request to the API.

#### kill(wait_timeout: float | int | None = 1, verbosity_level: int = 2) → bool | None

Request a task to be killed.

This method requests that the current task is remotely killed.
If wait_timeout is None (default), the kill request is sent to the
backend and the method returns. However, if wait_timeout is
a positive number, the method waits up to wait_timeout seconds
to ensure that the task transitions to the KILLED state.
:param wait_timeout: Optional - number of seconds to wait
:type wait_timeout: int, float
:param for the kill command or None if only the request is to be sent.:
:param verbosity_level: Optional. the verbosity of the logs when the
:type verbosity_level: int
:param task signal is sent and when the task is killed. Verbosity 0:
:param produces no outputs:
:type produces no outputs: Default
:param 1 produces minimal outputs:
:type 1 produces minimal outputs: Default
:param and 2:
:type and 2: Default
:param produces extensive outputs.:

* **Returns:**
  - None if wait_timeout is None and the kill request was
    successfully sent;
  - True if wait_timeout> 0 and the task successfully transitioned
    to the KILLED state within wait_timeout seconds;
  - False if wait_timeout > 0 but the task didn’t transition
    to the KILLED state within wait_timeout seconds;

#### last_modified_file()

Display the last modified file for a given task.

This function retrieves and prints information about the most recently
modified file associated with a specified task. It validates that the
task computation has started before proceeding. If the task is invalid
or not started, an error message is printed to stderr.

#### list_files() → Tuple[str | None, int]

List the files in the task’s working directory.

This method will list the files, in real time, in the task’s working
directory. It will also print the files in a tree-like structure.

* **Returns:**
  A string with the formatted directory listing.
  The return code for the command. 0 if successful, 1 if failed.

#### print_summary(fhandle=<_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>)

#### remove_remote_files(verbose: bool = True) → bool

Removes all files associated with the task from remote storage.

* **Returns:**
  True if the files were removed successfully, False otherwise.

#### set_metadata(metadata: Dict[str, str])

Set metadata for the task.

Metadata is stored as key-value pairs, where both
keys and values must be strings.
Metadata can be useful for categorizing, searching,
and filtering tasks.

Example usage:

> task = simulator.run(…)
> # Add experiment information to the task
> task.set_metadata({

> > “study”: “study_1”,
> > “experiment”: “experiment_1”,
> > “description”: “This is a test experiment”,
> > “parameters”: “param1=1,param2=2”,

> })
* **Parameters:**
  **metadata** – A dictionary with the metadata to set.

#### *property* summary *: str*

It returns cached information about the task summary.

#### sync_context()

Enter context manager for blocking sync execution.

This turns an asynchronous task into a blocking task.
If an exception/ctrl+c is caught while in the context manager, the
remote task is killed.

Usage:
: task = scenario.simulate(…)
  <br/>
  with task.sync_context():
  : # If an exception happens here or ctrl+c is pressed, the task
    # will be killed.
    task.wait()

#### tail_files(tail_files: ~typing.List[str], lines: int, follow: bool, fout: ~typing.TextIO = <_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>, wait: bool = False)

Prints the result of tailing a list of files.

* **Parameters:**
  * **tail_files** – A list of files to tail.
  * **lines** – The number of lines to print.
  * **follow** – Whether to keep tailing a file or not. If True, tail_files
    will keep printing the new lines in the selected files as they
    are changed in real time. If False, it will print the tail and
    end.
  * **fout** – The file object to print the result to. Default is stdout.
  * **wait** – If True, the method will wait for the files to be created
    before tailing them.

#### wait(polling_period: int = 1, silent_mode: bool = False, download_std_on_completion: bool = True) → TaskStatusCode

Wait for the task to complete.

This method issues requests to the API.

* **Parameters:**
  * **polling_period** – How often to poll the API for the task status.
  * **silent_mode** – If True, do not print the task logs (stdout and stderr)
    to the console.
  * **download_std_on_completion** – Request immediate download of the
    standard files (stdout and stderr) after the task completes.
* **Returns:**
  The final status of the task.

#### wait_for_status(status: str, polling_period: int = 1, silent_mode: bool = False) → TaskStatusCode

Wait for the task to reach a specific status or complete.

This method issues requests to the API.

* **Parameters:**
  * **polling_period** – How often to poll the API for the task status.
  * **silent_mode** – If True, do not print to stdout.
  * **status** – Return when the task reaches the set status or if the
    task reaches a terminal status.
* **Returns:**
  The final status of the task.

::docsbannersmall
::
