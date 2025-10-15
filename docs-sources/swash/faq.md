**Find answers to commonly asked questions about SWAN.**

<br>

# FAQ

## 1. My simulation runs almost to the end, but then fails and the status changes to “Failed.” Why?
One of the most common reasons for this issue is a mismatch between the version of SWASH you are requesting Inductiva to run (e.g., `11.01`) and the version for which your configuration file was originally created (e.g., `9.01`).

You can confirm this by checking your `stdout.txt` file. You might see that the simulation appears to be progressing normally, such as:

```
...
[ 91%]   step 23304
[ 91%]   step 23305
[ 91%]   step 23306
[ 91%]   step 23307
[ 91%]   step 23308
[ 91%]   step 23309
[ 91%]   step 23310
[ 91%]   step 23311
[ 91%]   step 23312
```

However, in your `stderr.txt` file, you’ll likely find an error similar to this:

```
--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 38 in communicator MPI_COMM_WORLD
with errorcode 4.

NOTE: invoking MPI_ABORT causes Open MPI to kill all MPI processes.
You may or may not see output from other processes, depending on
exactly when Open MPI kills them.
--------------------------------------------------------------------------
[api-6jn2hhpd76b7svd9itax53f6y-2twc:02531] 55 more processes have sent help message help-mpi-api.txt / mpi-abort
[api-6jn2hhpd76b7svd9itax53f6y-2twc:02531] Set MCA parameter "orte_base_help_aggregate" to 0 to see all help / error messages
```

This typically indicates that your configuration file was created for an older SWASH version.

**Solution**: Make sure you explicitly specify the correct version when creating the SWASH object. For example:

```python
# Initialize the Simulator
swash = inductiva.simulators.SWASH(version="10.05")
```

You can check the [list of currently available versions](versions-and-containers.md) to ensure compatibility.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
