# Storage and Data Flow

When talking about concepts like storage at Inductiva, there are 3 levels you 
need to consider:

**Local Storage:** this is the storage of your local machine, from where you are 
typically writing Python scripts that call the Inductiva API. Inductiva does not 
have direct access to this storage. However, the API has primitives that allow you 
to exchange information from your local storage to your Personal Remote Storage 
(see below).

**Personal Remote Storage:** this is a folder that lives on Inductiva cloud storage 
that is exclusively dedicated to you, and can only be accessed by you. We use this 
space both to store input files that you submit via the API (typically when you 
call the run() method of a simulator object), and to share with you the potentially 
large output files generated by the simulators you invoke.

**Worker Storage:** this is the storage that is available on the worker VMs that pick 
up simulation tasks. This is used to store the input files required for the simulation 
to run as well as the files produced by the simulation before these get transferred 
to your Personal Remote Storage. You don’t have direct access to this storage. All 
communication to your Local Storage needs to be done via your Personal Remote Storage folder.

So, what is the typical flow of that when you invoke a remote simulator using the 
Inductiva API? 

## Local Storage

Let’s start by assuming that somewhere in your local storage, typically in a 
folder dedicated to your project, you have several files that are required for 
running the simulation. Usually this data includes one or more files describing 
the simulation case, and file specifying how the simulator should be parameterized, 
as well as files describes assets (e.g. 3D shapes, information about chemical compounds, 
bathymetric profiles, etc) that will be used in or that are themselves the target 
of the simulation. All this data is input to the simulator software and, of course, 
will have to be sent to Inductiva machines where the simulator is going to execute. 
So, the first step in the data flow is to upload all these files to our server. 
This upload is triggered when you call the run method of the simulator object 
you are using. 

Here is an example. Let us assume you are developing a coastal dynamics study 
using Reef3D, and you have all the required input files and assets stored in the 
subdirectory `my_input_data_dir` located inside your project folder on your local 
machine. The following piece of code illustrates this situation:

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Initialize the simulator object
simulator = inductiva.simulators.REEF3D()

# Invoke the run() method of the simulator object. 
# This will trigger the packing and uploading the data
task = simulator.run(input_dir="my_input_data_dir",
                     on=machine_group)

# Terminate the machine group
machine_group.terminate()
```

The moment you invoke `run()` you start the uploading process. The folder `my_input_data_dir`
is zipped and the corresponding zip file is uploaded to Inductiva servers.
 
You can check what is happening when you invoke a simulator via the API. 
If you look at the logs produced at you will be able to see a message like this 
right in the beginning of the process execution:

```bash
Task Information:
> ID:                    tc7cwuer45kfzuw8t93r6dxa8
> Method:                swash
> Local input directory: swash-resources-example
> Submitting to the following computational resources:
 >> Default queue with c2-standard-4 machines.
Preparing upload of the local input directory swash-resources-example (160 B).
Local input directory successfully uploaded.
Task tc7cwuer45kfzuw8t93r6dxa8 submitted to the default queue.
```

## Personal Remote Storage

Once the zip file gets to the Inductiva server, it is immediately 
transferred to your Personal Remote Storage area, under a folder whose 
name is, by default, the ID for the simulation task you invoked. You can 
check the contents of your  Personal Remote Storage programmatically via 
the API or by using the CLI. Next, we show how you would be able to check the uploaded zip file using the CLI.

To check your personal storage area, you can do a general listing of the contents with:
```bash
$ inductiva storage ls

       NAME                             SIZE          CREATION TIME
       tc7cwuer45kfzuw8t93r6dxa8/       1.53 MB       08 Feb, 14:08:44
       qetbcydbymfg9r3eqri7jbekh/       7.36 MB       08 Feb, 14:07:34
       sk9zbcdkfo0124tw0jvro8if0/       7.36 MB       08 Feb, 14:07:33
       osd7r4onvetxmrxggpohc4wdc/       7.26 MB       08 Feb, 14:07:33
       gpawvd4qbnq36vhy3z0kddj7u/       7.28 MB       08 Feb, 14:07:32
       qtjs7n5xnaixfuhu8jm03xv38/       7.36 MB       08 Feb, 14:07:32
       9e3bgdgwqahrwaxvrk06q5mhn/       7.36 MB       08 Feb, 14:07:31
       5ry5h8q26o0fxqs8mymsy0d7r/       7.27 MB       08 Feb, 14:07:30
       jbxo7dc9pypqxzm53mw0jjk5p/       7.29 MB       08 Feb, 14:07:30
       9uhbxuzy2bqnjyyt4arxheqwc/       7.36 MB       08 Feb, 14:07:29
       i2ge334hdy4kinwvmau5dtwxx/       7.36 MB       08 Feb, 14:07:29

```

The simulation we have just invoked has the task ID `hzgk5ngzk28a39qa7mesv0snk`
and we can check that its contents were correctly submitted to the server by listing
the specific contents of the task folder with:

```bash
$ inductiva storage ls tc7cwuer45kfzuw8t93r6dxa8

       NAME             SIZE          CREATION TIME
       input.zip        1.06 KB       08 Feb, 14:08:44
                        0 B           08 Feb, 14:08:44
```

Once your simulation task gets picked up by a Worker, its input files need to be
transferred from your Personal Remote Storage to the corresponding VM. Typically,
this VM lives in the same region of the Google Cloud storage, and so moving data
is pretty fast.

## Worker Storage

Of course, the receiving VM needs to have enough storage space to execute your simulation. 
Typically, the input data for a simulation is relatively small. In the example above, the 
files required to run the simulation only have `1.53` MB. The 
challenge is the size of the outputs produced by running the simulation, which can easily 
get to dozens of GB, so VMs need to have large enough disk space installed.

Now, VM storage space can turn out to be pretty expensive, so we allow users to 
explicitly define the amount of VM storage dedicated to storing the results of the 
simulation, taking into account what they believe is the reasonable effective and
realistic amount needed. The size of the storage in the computational resources
is selected in the initialization with the parameter `data_disk_gb`. Note that when
running simulations on the shared pool of resources it is not possible to configure
the disk storage, therefore, simulations that may occupy more than the available,
30 GB, will fail.

When the simulation finishes running, the output files are uploaded to the respective
task folder in the user's remote storage and they become available to be downloaded
by the user whenever required. When this process finishes the corresponding data 
in the worker is freed. Therefore, the worker storage only needs to account for a
simulation at a time.
