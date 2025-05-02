# Run Your First Simulation
This tutorial will show you how to run COAWST simulations using the Inductiva API. 

## Compilation Requirements
Unlike many simulators, COAWST must be compiled for each specific configuration. This means that you cannot use a single pre-compiled version across simulations - COAWST must be compiled using the appropriate settings 
tailored to your use case.

To simplify this process, you'll need to include a few extra files in your simulation inputs, alongside your standard configuration and data files:
- A **COAWST build script** (typically named `build_coawst.sh`)
- The **header files** necessary for compilation
- Any additional files specific to your setup (e.g., custom switch files)

If you are using files already provided in the COAWST repository, you don’t need to upload them along with your input files. 
Instead, simply configure your build script to reference the path where Inductiva exposes those files. Alternatively, you can copy the 
necessary files from the exposed COAWST directory into your input files (more on this in the following section).

For each simulation, the COAWST directory will be available at: 
📂 `/workdir/output/artifacts/__COAWST`  

You are free to access and use any files within this directory.

Additionally, all input files will be located at:
📂 `/workdir/output/artifacts/`

Please keep this in mind when working with absolute paths.

### Initialization Commands  
When running COAWST, you can specify `init_commands` - a set of commands executed
before compilation. These are useful for copying the necessary files from the COAWST
directory to your working directory.  

For example:
```python
init_commands = [
    # Copy LANDUSE.TBL for the simulation to ".".
    # Note that "." points to /workdir/output/artifacts/
    "cp /workdir/output/artifacts/__COAWST/LANDUSE.TBL ."
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

## Run a COAWST Simulation With ROMS, SWAN and WRF

### Objective
We will cover the `JOE_TC/DiffGrid` use case from the [COAWST GitHub repository](https://github.com/DOI-USGS/COAWST), where
- **ROMS** and **SWAN** share the same grid
- **WRF** runs on a separate grid.

This configuration illustrates how to set up and execute models using different grid resolutions within COAWST.

### Prerequisites  
1. Download the **JOE_TC/DiffGrid** project from [this link](https://github.com/DOI-USGS/COAWST/tree/f1a4250bc64bf0c4f9d521effb47d85837c92e8a/Projects/JOE_TC/DiffGrid)
2. Download the standard **`build_coawst.sh`** script from [this link](https://github.com/DOI-USGS/COAWST/blob/f1a4250bc64bf0c4f9d521effb47d85837c92e8a)
3. Download the following files from the [JOE_TC folder](https://github.com/DOI-USGS/COAWST/tree/main/Projects/JOE_TC)
   - `namelist.input`
   - `wrfbdy_d01`
   - `wrfinput_d01`
 
Create a folder named `JOE_TC_DiffGrid` and place all files inside. Your folder structure should look like this:

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

### Update Your Input Files  
In this section, you will update the input files (`build_coawst.sh` and `.in` files) to ensure they reference the correct paths.

#### Update `build_coawst.sh`
Make the following changes:

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
These are all the necessary modifications to the script. Once updated, your `build_coawst.sh` will be properly 
configured for the compilation process.

#### Update the `.in` Files
Next, you need to update your simulation files to reflect the correct paths used in the COAWST repository. 
Typically, this involves modifying references from `Projects/JOE_TC/DiffGrid/file.txt` to just `file.txt`.

Start by updating the **coupling_joe_tc.in** file:
1. Change `WAV_name` to `WAV_name = INPUT_JOE_TC_COARSE`
2. Change `OCN_name` to `OCN_name = ocean_joe_tc_coarse.in`
3. Change `SCRIP_COAWST_NAME` to `SCRIP_COAWST_NAME = scrip_joe_tc_diffgrid.nc`
4. Change `W2ONAME` to `W2ONAME == wav1_to_ocn1_weights.nc`
5. Change `W2ANAME` to `W2ANAME == wav1_to_atm1_weights.nc`
6. Change `A2ONAME` to `A2ONAME == atm1_to_ocn1_weights.nc`
7. Change `A2WNAME` to `A2WNAME == atm1_to_wav1_weights.nc`
8. Change `O2ANAME` to `O2ANAME == ocn1_to_atm1_weights.nc`
9. Change `O2WNAME` to `O2WNAME == ocn1_to_wav1_weights.nc`

Next, update the wave model file `INPUT_JOE_TC_COARSE`:
1. Change `READGRID COORDINATES 1 'Projects/JOE_TC/DiffGrid/joe_tc_coarse_grid_coord.grd' 4 0 0 FREE`
to `READGRID COORDINATES 1 'joe_tc_coarse_grid_coord.grd' 4 0 0 FREE`
2. Change `READINP BOTTOM 1 'Projects/JOE_TC/DiffGrid/joe_tc_coarse_bathy.bot' 4 0 FREE `
to `READINP BOTTOM 1 'joe_tc_coarse_bathy.bot' 4 0 FREE `

Lastly, update the ocean model file `ocean_joe_tc_coarse.in`:
1. Change `VARNAME` to `VARNAME = /workdir/output/artifacts/__COAWST/ROMS/External/varinfo.dat`
2. Change `GRDNAME`to `GRDNAME == joe_tc_coarse_grd.nc`
3. Change `ININAME` to `ININAME == joe_tc_coarse_ocean_init.nc`

### Run Your Simulation  
You're now ready to send your simulation to the Cloud!

Here is the code required to run a COAWST simulation using the Inductiva API:

```python
"""COAWST Simulation."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2-standard-4",
    spot=True
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

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
task.print_summary()
```

**Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

In this example, we're using a relatively small cloud machine (`c2-standard-4`), which is equipped with 4 virtual CPUs. 
COAWST requires a precise core allocation for its simulations, meaning the number of CPUs must exactly match the simulation's configuration. 
In this case, the input files indicate that the simulation should run on **3 cores** (`n_vcpus=3`) - no more, no less.

This configuration is defined in the `coupling_joe_tc.in` file:

```
! Number of parallel nodes assigned to each model in the coupled system.
! Their sum must be equal to the total number of processors.

   NnodesATM =  1                    ! atmospheric model
   NnodesWAV =  1                    ! wave model
   NnodesOCN =  1                    ! ocean model
   NnodesHYD =  0                    ! hydrology model
```

Each component of the simulation is assigned a specific number of cores. 
While we can increase this number later, we'll keep it as is for now.

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
		├> 1.234 s         rm -r __COAWST
    └> 1.065 s         clean_all_sim_links
	Finalizing                at 28/02, 05:43:05      132.651 s
	Success                   at 28/02, 05:45:18      

Data:
	Size of zipped output:    4.45 GB
	Size of unzipped output:  4.79 GB
	Number of output files:   43

Estimated computation cost (US$): 0.71 US$
```

The simulation details might seem complex initially, but let's focus on the `In Progress` stage, as this is the part specific to your simulation. 
All other steps are common to every simulation run on Inductiva.

### Understanding the `In Progress` steps
The "In Progress" section outlines the commands executed during your simulation, along with their durations. Below is a breakdown of the key steps:

- **`cp -r /opt/COAWST /workdir/output/artifacts/__COAWST`**  
  - Copies the COAWST directory into the input files directory, allowing you to compile and access all files related to COAWST.

- **`create_all_sim_links`**  
  - Creates symbolic links (shortcuts) in the working directory for all COAWST files. This is necessary because some simulations require specific files
  (e.g., `CAMtr_volume_mixing_ratio`, `LANDUSE.TBL`) to be available 
  in the working directory. By creating these symbolic links, we avoid the need to send all those files with the input.
  - If a file with the same name already exists, the link is skipped, allowing you to provide your own version. For example, if you include your own
  `LANDUSE.TBL` file in the input directory, it will be used for the simulation instead of the default version from the COAWST folder.

- **`bash build_coawst.sh`**  
  - Runs your provided script to build COAWST.  
  - This script is executed inside `/workdir/output/artifacts`, and COAWST is
  located at `/workdir/output/artifacts/__COAWST`.  

- **`coawstM coupling_joe_tc.in`**  
  - Runs the simulation in parallel.

- **`rm -r __COAWST`**  
  - Cleans up the simulation directory by removing the COAWST folder, reducing the overall output size.

- **`clean_all_sim_links`**
  - Cleans all the symbolic links created by the `create_all_sim_links` command.

These steps ensure that your COAWST simulation runs efficiently and is well-managed in the cloud. 🚀

### Scale Up Your Simulation
Based on the execution times, the compilation took
**1,335 seconds** (around **22 minutes**), while the simulation itself ran for
**35,278 seconds** (approximately **9 hours and 47 minutes**).  

It's important to note that compilation time varies depending on the
**COAWST configuration** chosen. However, the simulation runtime is
primarily influenced by the **computational resources allocated**—including the
number of virtual CPUs.  

In this section, we'll explore strategies to **scale up your simulation**,
in order to reduce the simulation time.

#### Update Your Input Files
As mentioned earlier, the number of virtual CPUs used for your simulation must exactly match the configuration specified in the input files. 
Therefore, to scale up the simulation, you'll need to modify the following three files:

- `coupling_joe_tc.in`
  - Modify `NnodesATM`: Number of virtual CPUs assigned to the atmospheric model.
  - Modify `NnodesWAV`: Number of virtual CPUs assigned to the wave model.
  - Modify `NnodesOCN`: Number of virtual CPUs assigned to the ocean model.
  - **Restriction**: The sum of `NnodesATM + NnodesWAV + NnodesOCN` must equal the `n_vcpus` passed
  in the python script.
- `ocean_joe_tc_coarse.in`
  - Modify `NtileI`: I-direction partition.
  - Modify `NtileJ`: J-direction partition.
  - **Restriction**: The product of `NtileI * NtileJ` must equal `NnodesOCN` (as defined in the `coupling` file).
- `namelist.input`
  - `nproc_x`
  - `nproc_y`
  - Set `nproc_x * nproc_y` equal to `NnodesATM` to take full
  advantage of the virtual cores assigned to the atmosferic model.

Below is a list of the simulations with their respective results:

|   Machine Type  | Virtual CPUs |     Execution Time     |   Cost   |
|:---------------:|:------------:|:----------------------:|:--------:|
|  c2-standard-4  |       4      | 9 hours and 47 minutes | 0.71 US$ |
|  c2-standard-60 |      60      |  1 hour and 3 seconds  | 1.37 US$ |

For these simulations, we divided the number of virtual CPUs available on each machine equally among the three models (1, 20, and 36 virtual CPUs, respectively).
We used the following values for `Ntile` and `nproc`: (1 1), (4 5), and (6 6), respectively.

## Conclusion
We’ve covered the key steps for setting up and running a **COAWST** simulation using the **Inductiva API**. We also explored the necessary modifications
to input files for compiling and running **COAWST**.

By following this guide, you now have a clearer understanding of how to configure and efficiently run COAWST simulations on Inductiva's platform.

Happy simulations! 🚀
