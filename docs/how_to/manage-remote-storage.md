<meta http-equiv="refresh" content="0; url=https://tutorials.inductiva.ai/how_to/manage-remote-storage.html">
<script type="text/javascript">
    window.location.href = "https://tutorials.inductiva.ai/how_to/manage-remote-storage.html";
</script>

# Manage Remote Storage

The Inductiva API provides you with tools to manage your remote storage effectively. 
With the Inductiva storage module, you can easily navigate your storage, evaluate
the space used, and delete specific directories as needed.
Let's illustrate this with some examples.

## Determining the amount of storage in use

As you start running simulations, your remote personal storage will get filled up.
Hence, at times it will be useful to monitor the amount of storage space being utilized.
This can be achieved as follows:

### Python
```python
import inductiva
space_used = inductiva.storage.get_space_used()
```
### CLI
```bash
$ inductiva storage size
```

In both cases, we receive the following message:

```
Total user's remote storage in use: 32.4 GB
```

## Viewing storage contents

After determining the total storage space used, you may want to identify which directories are consuming the most storage. To be more specific, you may list
the contents of your directory and sort them by size or creation date as follows:

**Python**

```python
import inductiva
inductiva.storage.listdir(max_results=10, order_by="size", sort_order="desc")
```

**CLI**

```bash
$ inductiva storage ls --max-results 10 --order-by size --sort-order desc
```

In this way, we obtain a listing of the 10 largest directories within our user's
remote storage.
```bash

       NAME                             SIZE           CREATION TIME
       0bet8jrpp2gz974n42nsd9n2p/       56.11 MB       06 Feb, 11:32:29
       f8joznwc9xf9a4nypcaei6v2s/       12.79 MB       07 Feb, 09:16:55
       r4kerxf4b53krgn0s3fyece3b/       11.92 MB       07 Feb, 11:47:48
       j9qzrpiohgt7x97od3tw4wccd/       11.74 MB       07 Feb, 11:47:46
       iqi71gonoacfj7fknox3rvnq2/       11.52 MB       07 Feb, 11:47:45
       dxmnxdrfrv84pfbzbvm9v0dat/       11.43 MB       07 Feb, 11:47:43
       bgtwgnnyq5qa5hecegzdx6okr/       11.36 MB       07 Feb, 11:47:40
       6a2h1wnxywpea8jfoxiikdjf7/       1.53 MB        07 Feb, 13:47:03
       hkrr6qtiuu8uatwgrydc1e6q5/       1.53 MB        07 Feb, 09:26:10
       97ujsu1uc48oc3xi1572dj1uh/       1.53 MB        07 Feb, 11:31:30
```

In case, you want to be more specific and examine the contents of a specific folder,
you can pass a path to the `listdir` method and/or the CLI subcommand as follows:

**Python**

```python
import inductiva
inductiva.storage.listdir(path = "0bet8jrpp2gz974n42nsd9n2p")
```

**CLI**

```bash
$ inductiva storage ls 0bet8jrpp2gz974n42nsd9n2p
```

Outputs the following::
```bash
       NAME             SIZE           CREATION TIME
       output.zip       52.37 MB       06 Feb, 11:34:26
       input.zip        3.74 MB        06 Feb, 11:32:32
                        0 B            06 Feb, 11:32:29
```

### Removing directories

Whenever space needs to be freed up, or in general, when you want to remove a directory
from your remote storage, you can use the task interface to do it. 

The table above provides valuable information that can guide your decision to remove certain directories. 

**Python**

```python
import inductiva
inductiva.tasks.Task(task_id).remove_remote_files()
```

**CLI**

```console
$ inductiva storage rm reef3d_simulation
You are about to remove the following paths from your remote storage space:
  - reef3d_simulation
Are you sure you want to proceed (y/[N])? y
Removing reef3d_simulation in the user's remote storage.
Successfully removed remote path 'reef3d_simulation'.
```

Notice that you get prompted, to make sure you want to proceed with the removal.
In case you wish to remove everything, use the flag `--all`, like `inductiva storage rm --all`.

In the end, we can verify that the folder we have named `reef3d_simulation` is removed from the user's remote storage and that nothing appears in the listing.

```bash
$ inductiva storage ls reef3d_simulation

       NAME         SIZE         CREATION TIME

```
