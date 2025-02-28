This tutorial will guide you through the process of setting up and running a
**COAWST** simulation using the **Inductiva API**.

# COAWST

**COAWST (Coupled Ocean-Atmosphere-Wave-Sediment Transport)** modeling system is
a powerful numerical simulation framework designed to study coastal and oceanic
processes. It integrates multiple modeling components, including:

- **ROMS (Regional Ocean Modeling System)** for ocean dynamics  
- **WRF (Weather Research and Forecasting Model)** for atmospheric processes  
- **SWAN (Simulating Waves Nearshore)** for wave dynamics  
- A **sediment transport module** for coastal sediment simulations  

By coupling these models, COAWST allows researchers to analyze complex
interactions between the ocean, atmosphere, and coastal environments. It is
widely used for applications such as
**storm surge prediction, sediment transport studies, and climate impact assessments**.

## Supported Versions

The Inductiva API currently supports the following COAWST version:  
- **v3.8**

## Compilation Requirements

Unlike many other simulators, COAWST requires a compilation specifically for
each configuration. This means that a single compiled version of
COAWST cannot be used for all simulations—you must compile it with the
appropriate settings for your use case.

To simplify this process, we require a few additional files as part of your
simulation input, in addition to the usual configuration and data files:

- A **COAWST build script**, typically named `build_coawst.sh` or something similar.
- Any **header files** necessary for compiling COAWST.
- Any other file that you need for your simulation.

If you are using standard switch and header files provided in the COAWST
repository, you don't need to include them manually. Instead, configure your
build script to point to the correct paths.

For each simulation, the COAWST directory will be available at:  
📂 `/workdir/output/artifacts/__COAWST`  

You can freely access and use files within this directory. Additionally, all
input files will be placed in:  
📂 `/workdir/output/artifacts/`

Keep this in mind when using absolute paths.  

### Initialization Commands  

When running COAWST, you can specify `init_commands`, a set of commands executed
before compilation. These are useful for copying necessary files from the COAWST
directory to your working directory.  

#### Example Usage  

```python
init_commands = [
    # Copy LANDUSE.TBL for the simulation
    "cp /workdir/output/artifacts/__COAWST/LANDUSE.TBL ."
    # Or (same result as the last command)
    "cp /workdir/output/artifacts/__COAWST/LANDUSE.TBL /workdir/output/artifacts/LANDUSE.TBL"
]

# Run simulation
task = coawst.run(
    input_dir="/Path/to/input_files",
    sim_config_filename="sim_config.in",
    build_coawst_script="build_coawst.sh",
    init_commands=init_commands,  # Pre-compilation commands
    n_vcpus=360,
    use_hwthread=True,
    on=cloud_machine
)
```

Follow the next sections to run a COAWST simulation with Inductiva.

## Running a COAWST Simulation With ROMS, SWAN and WRF

### Objective

In this simulation, we will run the **JOE_TC/DiffGrid** case using **COAWST**,
where **ROMS** and **SWAN** share the same grid, while **WRF** runs on a
different grid. This setup demonstrates how to configure and execute models with
different grid resolutions within COAWST.

The goal is to show you how you need to configure your build script in order to run
your simulation with us.


### Prerequisites  

To get started, download the **JOE_TC/DiffGrid** project from the
[Official COAWST GitHub Repository](https://github.com/DOI-USGS/COAWST/tree/f1a4250bc64bf0c4f9d521effb47d85837c92e8a/Projects/JOE_TC/DiffGrid).  

You will also need the standard **`build_coawst.sh`** script which you can
download from [this link](https://github.com/DOI-USGS/COAWST/blob/f1a4250bc64bf0c4f9d521effb47d85837c92e8a).

Lastly, download the three files (`namelist.input`, `wrfbdy_d01` and `wrfinput_d01`)
present [in the JOE_TC folder](https://github.com/DOI-USGS/COAWST/tree/main/Projects/JOE_TC).

Once downloaded, place all files inside a folder named **`JOE_TC_DiffGrid`**.  

Your folder structure should now look like this:  
```
-rwxr-xr-x@ 1 paulobarbosa  staff      3085 Feb 27 10:47 INPUT_JOE_TC_COARSE
-rwxr-xr-x@ 1 paulobarbosa  staff     15528 Feb 27 10:56 build_coawst.sh
-rwxr-xr-x@ 1 paulobarbosa  staff      7266 Feb 28 07:12 coupling_joe_tc.in
-rwxr-xr-x@ 1 paulobarbosa  staff      2978 Nov 11 20:45 joe_tc.h
-rwxr-xr-x@ 1 paulobarbosa  staff    116775 Nov 11 20:45 joe_tc_coarse_bathy.bot
-rwxr-xr-x@ 1 paulobarbosa  staff   1671068 Nov 11 20:45 joe_tc_coarse_grd.nc
-rwxr-xr-x@ 1 paulobarbosa  staff    195000 Nov 11 20:45 joe_tc_coarse_grid_coord.grd
-rwxr-xr-x@ 1 paulobarbosa  staff   4763160 Nov 11 20:45 joe_tc_coarse_ocean_init.nc
-rwxr-xr-x@ 1 paulobarbosa  staff      5132 Feb 28 07:19 namelist.input
-rwxr-xr-x@ 1 paulobarbosa  staff    167380 Feb 28 07:13 ocean_joe_tc_coarse.in
-rwxr-xr-x@ 1 paulobarbosa  staff  25403104 Nov 11 20:45 scrip_joe_tc_diffgrid.nc
-rwxr-xr-x@ 1 paulobarbosa  staff  46541404 Nov 11 20:45 wrfbdy_d01
-rwxr-xr-x@ 1 paulobarbosa  staff  70658632 Nov 11 20:45 wrfinput_d01
```

---

### Updating Your Input Files  

In this section, we will update the input files (`build_coawst.sh` and `.in`
files) to ensure they point to the correct paths.  

#### Updating the Build Script  

We need to modify a few settings in **`build_coawst.sh`** to ensure a proper setup.  

Here are the required changes:  

1. **Set the correct application name:**  
   - Update `COAWST_APPLICATION` to match your header file (`joe_tc.h`),
   capitalized and without the file extension:  
     ```bash
     export   COAWST_APPLICATION=JOE_TC
     ```  

2. **Set the root directory:**  
   - Update `MY_ROOT_DIR` to:  
     ```bash
     export   MY_ROOT_DIR=/workdir/output/artifacts/__COAWST
     ```  

3. **Specify the MPI implementation:**  
   - Change `which_MPI` to `openmpi`:  
     ```bash
     export   which_MPI=openmpi
     ```  

4. **Update header and analytical directories:**  
   - Ensure `MY_HEADER_DIR` and `MY_ANALYTICAL_DIR` point to the correct location
   where `joe_tc.h` is stored:  
     ```bash
     export   MY_HEADER_DIR=/workdir/output/artifacts
     export   MY_ANALYTICAL_DIR=/workdir/output/artifacts
     ```  

These are all the required modifications to the script. Once updated, your
**`build_coawst.sh`** will be correctly configured for the compilation process.

#### Updating the `.in` Files

Next, we need to update our simulation files to reflect the correct paths used
in the COAWST repository. Typically, this involves modifying any references
from `Projects/JOE_TC/DiffGrid/file.txt` to simply `file.txt`.

Starting with **coupling_joe_tc.in**:
1. Change `WAV_name` to `WAV_name = INPUT_JOE_TC_COARSE`
2. Change `OCN_name` to `OCN_name = ocean_joe_tc_coarse.in`
3. Change `SCRIP_COAWST_NAME` to `SCRIP_COAWST_NAME = scrip_joe_tc_diffgrid.nc`
4. Change `W2ONAME` to `W2ONAME == wav1_to_ocn1_weights.nc`
5. Change `W2ANAME` to `W2ANAME == wav1_to_atm1_weights.nc`
6. Change `A2ONAME` to `A2ONAME == atm1_to_ocn1_weights.nc`
7. Change `A2WNAME` to `A2WNAME == atm1_to_wav1_weights.nc`
8. Change `O2ANAME` to `O2ANAME == ocn1_to_atm1_weights.nc`
9. Change `O2WNAME` to `O2WNAME == ocn1_to_wav1_weights.nc`

We can now update the wave model file `INPUT_JOE_TC_COARSE`:
1. Change `READGRID COORDINATES 1 'Projects/JOE_TC/DiffGrid/joe_tc_coarse_grid_coord.grd' 4 0 0 FREE`
to `READGRID COORDINATES 1 'joe_tc_coarse_grid_coord.grd' 4 0 0 FREE`
2. Change `READINP BOTTOM 1 'Projects/JOE_TC/DiffGrid/joe_tc_coarse_bathy.bot' 4 0 FREE `
to `READINP BOTTOM 1 'joe_tc_coarse_bathy.bot' 4 0 FREE `

Lastly, we can update the ocean model file `ocean_joe_tc_coarse.in`:
1. Change `VARNAME` to `VARNAME = /workdir/output/artifacts/__COAWST/ROMS/External/varinfo.dat`
2. Change `GRDNAME`to `GRDNAME == joe_tc_coarse_grd.nc`
3. Change `ININAME` to `ININAME == joe_tc_coarse_ocean_init.nc`


### Running Your Simulation  

Now that you've made all the necessary changes to the input files, it's time to
run your simulation. Below is the Python code you need to execute:  

```python
"""COAWST Simulation."""
import inductiva

# Instantiate cloud machine
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2-standard-4"
)

# Initialize the simulator
coawst = inductiva.simulators.COAWST(
    version="3.8")

# Run simulation
task = coawst.run(
    input_dir="/Path/to/JOE_TC_DiffGrid",
    sim_config_filename="coupling_joe_tc.in",
    build_coawst_script="build_coawst.sh",
    n_vcpus=3,
    on=cloud_machine
)

task.wait()
cloud_machine.terminate()

task.download_outputs()
task.print_summary()
```

In this example, we're using a relatively small cloud machine (`c2-standard-4`),
which has 4 virtual CPUs. COAWST requires a strict core allocation for its
simulations, meaning the number of CPUs must exactly match the simulation's
configuration. Here, the input files specify that the simulation must run on
**3 cores**, no more, no less (`n_vcpus=3`).  

This is defined in the `coupling_joe_tc.in` file:  

```
! Number of parallel nodes assigned to each model in the coupled system.
! Their sum must be equal to the total number of processors.

   NnodesATM =  1                    ! atmospheric model
   NnodesWAV =  1                    ! wave model
   NnodesOCN =  1                    ! ocean model
   NnodesHYD =  0                    ! hydrology model
```

Each component of the simulation is assigned a specific number of cores. We can
increase this number later, but for now, we'll keep it as is.  

Once the simulation is complete, we terminate the machine, download the results
and print a summary of the simulation as shown below.

```
inductiva tasks info 6wt3dp49uhy45y708x848eu2y
Task status: Success

Timeline:
	Waiting for Input         at 27/02, 19:32:00      7.963 s
	In Queue                  at 27/02, 19:32:08      15.545 s
	Preparing to Compute      at 27/02, 19:32:24      11.374 s
	In Progress               at 27/02, 19:32:35      36630.294 s
		├> 12.082 s        cp -r /opt/COAWST /workdir/output/artifacts/__COAWST
		├> 1.063 s         create_all_sim_links
		├> 1335.318 s      bash build_coawst.sh
		├> 35278.638 s     /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 3 coawstM coupling_joe_tc.in
		└> 1.234 s         rm -r __COAWST
	Finalizing                at 28/02, 05:43:05      132.651 s
	Success                   at 28/02, 05:45:18      

Data:
	Size of zipped output:    4.45 GB
	Size of unzipped output:  4.79 GB
	Number of output files:   43

Estimated computation cost (US$): 0.71 US$

Go to https://console.inductiva.ai/tasks/6wt3dp49uhy45y708x848eu2y for more details.
```

The simulation details may seem overwhelming at first, but let's focus on the
`In Progress` stage, as this is the part unique to your simulation. All other
steps follow the same process for any simulation run on Inductiva.

### Understanding the `In Progress` steps

The **"In Progress"** section details the commands executed during your
simulation, along with their durations. Below is a breakdown of the key steps:

- **`cp -r /opt/COAWST /workdir/output/artifacts/__COAWST`**  
  - Copies the COAWST directory into the input files directory for you to compile
  and have access to all files related to COAWST.

- **`create_all_sim_links`**  
  - Creates symbolic links (shortcuts to other files) in the working directory
  for all COAWST files. We do this because some simulations require specific files
  (e.g., `CAMtr_volume_mixing_ratio`, `LANDUSE.TBL`) to be present in the working
  directory. By creating this symbolic links we avoid having to send all those files
  in the input files.
  - If a file with the same name already exists, the link is skipped, allowing
  you to provide your own version.  

- **`bash build_coawst.sh`**  
  - Runs your provided script to build COAWST.  
  - This script is executed inside `/workdir/output/artifacts`, and COAWST is
  located in `/workdir/output/artifacts/__COAWST`.  

- **`coawstM coupling_joe_tc.in`**  
  - Runs the simulation using `mpirun`, which automatically sets the number of
  cores based on the `n_vcpus` value used on the python script.

- **`rm -r __COAWST`**  
  - Cleans up the simulation directory by removing the COAWST folder,
  keeping the output size minimal.  

With these steps, your COAWST simulation is successfully executed and managed in the cloud. 🚀

## Scaling Up Your Simulation

Looking at the execution times, we can see that the compilation took
**1,335 seconds** (around **22 minutes**), while the simulation itself ran for
**35,278 seconds** (approximately **9 hours and 47 minutes**).  

It's important to note that compilation time varies depending on the
**COAWST configuration** you choose. However, the simulation runtime is
primarily influenced by the **computational resources allocated**—including the
number of virtual CPUs.  

In this section, we'll explore strategies to **scale up your simulation**,
in order to reduce the simulation time.


### Updating your input files

As we stated before, the number of virtual CPUs you use for your simulation needs
to match exctly with the configuration of the input files. So, in order to scale
up the simulation we need to edit 3 files:

- `coupling_joe_tc.in`
  - `NnodesATM`: number of virtual CPUs to assinged to the atmospheric model.
  - `NnodesWAV`: number of virtual CPUs to assinged to the wave model.
  - `NnodesOCN`: number of virtual CPUs to assinged to the ocean model.
  - `NnodesATM + NnodesWAV + NnodesOCN` needs to be equal to the `n_vcpus` passed
  in the python script.
- `ocean_joe_tc_coarse.in`
  - `NtileI`: I-direction partition.
  - `NtileJ`: J-direction partition.
  - `NtileI * NtileJ` needs to be equal to `NnodesOCN` (defined in the `coupling` file).
- `namelist.input`
  - `nproc_x`
  - `nproc_y`
  - We will make `nproc_x * nproc_y` equal to `NnodesATM` in order to take full
  advantage of the virtual cores assigned to the atmosferic model.

  Here is a small list of simulation with the respective configurations and the results:

  |  Machine Type  | Virtual CPUs |   NnodesATM  |  NnodesWAV | NnodesOCN | NtileI | NtileJ | nproc_x | nproc_y |     Execution Time     |   Cost   |
|:--------------:|:------------:|:------------:|:----------:|:---------:|:------:|:------:|:-------:|:-------:|:----------------------:|:--------:|
|  c2-standard-4 |       4      |       1      |      1     |     1     |    1   |    1   |    1    |    1    | 9 hours and 47 minutes | 0.71 US$ |
| c2-standard-60 |      60      |      20      |     20     |     20    |    4   |    5   |    4    |    5    |  1 hour and 3 seconds  | 1.37 US$ |
| c2-standard-16 |      16      | 10.2 seconds | 0.0011 US$ |           |        |        |         |         |                        |          |


In this tutorial, we covered the essential steps for setting up and running a
**COAWST** simulation using the **Inductiva API**. We explored the necessary
modifications to input files in order to compile and run **COAWST**.

By following this guide, you now have a better understanding of how to configure
and run COAWST simulations efficiently on Inductiva's platform.