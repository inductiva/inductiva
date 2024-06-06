# XBeach

XBeach is a simulator with a two-dimensional model for wave propagation,
sediment transport and morphological changes in the nearshore area. The
simulator is configured with a `params.txt` file that contains grid and
bathymetry info, wave input, flow input, morphological input, etc. in the form
of keyword/value pairs. If a `params.txt` cannot be found then XBeach will not
run. Other files are used to configure the grid and bathymetry profile, like
`bed.dep` for example, and other files with extra information that can be used
inside the `params.txt` to configure the simulator further.

We advise to always set the `mpiboundary` argument in the `params.txt` file, 
since we handle automatically the parallelization of the simulation, based on
the number of cores available in the machine.

## Example

```python
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "xbeach-input-example.zip", unzip=True)

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach()

# Run simulation with config files in the input directory
task = xbeach.run(input_dir=input_dir,
                  sim_config_filename="params.txt")

task.wait()
task.download_outputs()
```

## A more advanced example
We are now going to run a much longer simulation, namely one whose configuration 
scripts and input data are mae available via the
[GRIIDC](https://www.griidc.org/), a data repository based out of the Harte
Research Institute for Gulf of Mexico Studies at Texas A&M University-Corpus 
Christi. 

More specifically, we will run one of simulations described in 
"XBeach model setup and results for beach and dune enhancement scenarios on 
Galveston Island, Texas", which is made available
[here](https://data.griidc.org/data/HI.x833.000:0001).

We will start by downloading the required input files. At the bottom of the 
[dataset page](https://data.griidc.org/data/HI.x833.000:0001#individual-files)
you will find a widget for downloading indvidual files. We will need to 
download the contents of an entire directory. Using the widget, please
select: **Files >> XBeach_Model_Runs >> Beach_Nourish_Only >> Input**.

Select all files inside the directory using the checkbox and click on the 
elipsis "⋮" to download them:

![Select each of the files inside the Input folder](xbeach_source_files.png)

Move all files to a directory named  "Beach_Nourish_Only" inside your Inductiva
python project directory. Inside the  "Beach_Nourish_Only" you should then have
something like:

```
ls -las Beach_Nourish_Only 
total 130976
    0 drwxr-xr-x  12        384  5 Jun 16:13 .
    0 drwxr-xr-x   3         96  5 Jun 16:14 ..
    8 -rw-r--r--@  1       3069  5 Jun 16:12 README.txt
26184 -rw-r--r--@  1   13404906  5 Jun 16:12 bed.dep
26184 -rw-r--r--@  1   13404906  5 Jun 16:12 bedfricfile.txt
   16 -rw-r--r--@  1       5324  5 Jun 16:12 jonswap3.txt
26184 -rw-r--r--@  1   13404906  5 Jun 16:12 nebed.dep
    8 -rw-r--r--@  1       2281  5 Jun 16:12 params.txt
   16 -rw-r--r--@  1       4850  5 Jun 16:12 tide.txt
26184 -rw-r--r--@  1   13404906  5 Jun 16:12 x.grd
    8 -rw-r--r--@  1        635  5 Jun 16:12 xbeach.slurm
26184 -rw-r--r--@  1   13404906  5 Jun 16:12 y.grd
```

We are going to run this simulation, almost as it is, with only two small 
changes to the parametrization defined in ```params.txt```. First, we will need
to add an extra configuration parameter to account for the fact the API is 
already running XBeach v10+, but the original simulation configuration files
where prepared for an older version of XBeach. So will need to add the following
line to ```params.txt``` (you can add it right at the top, after the header):
```
single_dir = 0
``` 
Second, to reduce the time to complete this simulation, we will change the
stop time of simulation, ```tstop```, to just 34560 (10x less than original
parametrization). Find the ```tstop``` parameter in ```params.txt``` and change
the value to 34560.

That's it! We won't do any more changes. Let's start the simulation. The python
script we are going to use to trigger the simulation is:

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
    machine_type="c3d-highcpu-90",
    spot=True,
    data_disk_gb=20)
machine_group.start()

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach()

# Run simulation with config files in the input directory
task = xbeach.run(
    input_dir="Beach_Nourish_Only",
    sim_config_filename="params.txt",
    n_vcpus=90,
    on=machine_group)

task.wait()
task.print_summary()

# Terminate your dedicated MachineGroup at then end of the simulation.
machine_group.terminate()
```

There are a few of things you should notice in the script above:

* We are not requesting the outputs of the simulation to be donwloaded to our
local machine yet. That would be done by calling ```task.download_outputs()```
in this script right after ```task.wait()```. Instead, we opted for letting the
simulation finish and turning off the machine we are using. We will download the
outputs at a late time, after checking the size of the output that will be
automatically placed in our remote personal storage when the simulation ends.
* Finally, we are calling ```task.print_summary()``` that shows the times
spent at all stages of the process, including all auxiliary tasks, such as
moving data around between your local computer, your personal remote storage
space and the executer machine (i.e. the ```c3d-highcpu-90``` VM.)

Running this script, takes about 1h15, and should produce an output such as:

```
(base) lsarmento@Luiss-MacBook-Air xbeach % python run.py 
Registering MachineGroup configurations:
> Name:         api-6smqv1nofow4axk6vf99j1hwm
> Machine Type: c3d-highcpu-90
> Data disk size:    20 GB
> Number of machines: 1
> Spot:               True
> Estimated cloud cost of machine group: 1.230 $/h
Starting MachineGroup(name="api-6smqv1nofow4axk6vf99j1hwm"). This may take a few minutes.
Note that stopping this local process will not interrupt the creation of the machine group. Please wait...
Machine Group api-6smqv1nofow4axk6vf99j1hwm with c3d-highcpu-90 machines successfully started in 0:00:29.
Task Information:
> ID:                    dz3nbyekv0ds17hzi6227a7b9
> Method:                xbeach
> Local input directory: Beach_Nourish_Only
> Submitting to the following computational resources:
 >> Machine Group api-6smqv1nofow4axk6vf99j1hwm with c3d-highcpu-90 machines
Preparing upload of the local input directory Beach_Nourish_Only (67.04 MB).
Input archive size: 9.32 MB
Uploading input archive...
100%|██████████████████████████████████████████████████████████████████████████████| 9.32M/9.32M [00:10<00:00, 863kB/s]
Local input directory successfully uploaded.
Task dz3nbyekv0ds17hzi6227a7b9 submitted to the queue of the Machine Group api-6smqv1nofow4axk6vf99j1hwm with c3d-highcpu-90 machines.
Simulation metadata logged to: inductiva_output/task_metadata.json
Task dz3nbyekv0ds17hzi6227a7b9 configurations metadata saved to the tasks metadata file task_metadata.json in the current working directory.
Consider tracking the status of the task via CLI:
	inductiva tasks list --task-id dz3nbyekv0ds17hzi6227a7b9
Or, tracking the logs of the task via CLI:
	inductiva logs dz3nbyekv0ds17hzi6227a7b9
Task dz3nbyekv0ds17hzi6227a7b9 successfully queued and waiting to be picked-up for execution...
Task dz3nbyekv0ds17hzi6227a7b9 has started and is now running remotely.
Task dz3nbyekv0ds17hzi6227a7b9 completed successfully.
Wall time                         3263.68s
	input_upload:                   12.59
	container_image_download         1.02
	queue_time:                      8.53
	input_download:                  0.21
	input_decompression:             0.20
	computation:                  3217.49
	output_compression:             18.35
	output_upload:                   6.28
Data
	output_size:	398.82 MB
Terminating MachineGroup(name="api-6smqv1nofow4axk6vf99j1hwm"). This may take a few minutes.
Machine Group api-6smqv1nofow4axk6vf99j1hwm with c3d-highcpu-90 machines successfully terminated in 0:01:22.
```

The summary is pretty handy to understand that almost 99% of the (wall) time is
spent where is should be: on the computation stage, i.e. actually executing the
simulation.

Note: As seen in the code above we are using a machine with 90 vCPUs and,
in the method ```run()```, we are requesting the simulation to be parallelized
over all of those 90 vCPU. In some cases, parallelizing the simulation over only
half of the available vCPUs leads to better peformance. This is because the 
virtualization scheme of these VMs assigns two vCPU per underlying physical core
and so by setting ```n_vcpus``` to half the number of vCPUs we are implicitly
assinging one thread per physical core, which is many cases is more efficient.
However, this is NOT the case for this specific simulation with XBeach. In fact,
running the simulation on the same machine and setting ```n_vcpus =45``` will
make the computation about 35% slower. 

### Downloading simulation dats

Now, it is time to fecth the results. We will be downloading a zip file with
398.82MB of data (as shown in the summary). This can be done very conveniently
using the CLI. So, from your command line run (with the appropriate task id that
you can see above):

```
inductiva tasks download dz3nbyekv0ds17hzi6227a7b9
```

Depending on the speed of your internet connection, donwloading thie files may
take a few seconds or a few minutes:

```
Downloading simulation outputs to inductiva_output/dz3nbyekv0ds17hzi6227a7b9/output.zip.
100%|███████████████████████████████████████████████████████████████████████████████| 399M/399M [02:34<00:00, 2.59MB/s]
Uncompressing the outputs to inductiva_output/dz3nbyekv0ds17hzi6227a7b9.
```

If you now look under ```inductiva_output/dz3nbyekv0ds17hzi6227a7b9``` 
(please check the id of the task that you actually run on your terminal) inside
your project directory you should see something like:

```
ls -las inductiva_output/dz3nbyekv0ds17hzi6227a7b9
total 1313072
     0 drwxr-xr-x  29         928  6 Jun 07:09 .
     0 drwxr-xr-x   8         256  6 Jun 07:07 ..
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00001.bcf
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00002.bcf
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00003.bcf
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00004.bcf
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00005.bcf
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00006.bcf
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00007.bcf
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00008.bcf
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00009.bcf
 35224 -rw-r--r--   1    18033840  6 Jun 07:09 E_series00010.bcf
   256 -rw-r--r--   1      128821  6 Jun 07:09 XBlog.txt
     8 -rw-r--r--   1        1105  6 Jun 07:09 XBwarning.txt
     8 -rw-r--r--   1         860  6 Jun 07:09 ebcflist.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00001.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00002.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00003.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00004.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00005.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00006.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00007.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00008.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00009.bcf
 14096 -rw-r--r--   1     7213536  6 Jun 07:09 q_series00010.bcf
     8 -rw-r--r--   1         860  6 Jun 07:09 qbcflist.bcf
     8 -rw-r--r--   1        1248  6 Jun 07:09 stderr.txt
   256 -rw-r--r--   1      129208  6 Jun 07:09 stdout.txt
819328 -rw-r--r--   1   415526280  6 Jun 07:09 xboutput.nc
```

that contains all the files produced by XBeach, as well as two additional log 
files that the Inductiva API always captures: ```stdout.txt``` and 
```stderr.txt```. Now that you have your files locally, you can execute all 
sorts of post-processing steps on them, as you would if you had run your
simulations locally. 

That's it! You can now go bigger! You can start by trying to run the complete
simulation (the original parameter is ```tstop = 345600```) on an even faster 
machine such as a ```c3d-highcpu-360```!

Good luck!

## What to read next

If you are interested in XBeach, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [Reef3D](Reef3D.md)
* [SCHISM](SCHISM.md)
* [SWASH](SWASH.md)
* [SWAN](SWAN.md)
