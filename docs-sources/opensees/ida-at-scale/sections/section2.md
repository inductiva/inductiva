## Run 50 Simulations in Parallel
Let's walk through each of the six sections of the simulation script.

### Section 1: Allocating the Cloud Machine Group
We begin by allocating an Elastic Machine Group to run the simulations. This group is configured with 
a minimum and maximum number of active machines. Machines are dynamically turned on or off
as demanded, ensuring efficient resource utilization.

Since some simulations may take longer than others, this elastic setup prevents idle waiting. When 
a machine finishes its assigned tasks and no further simulations are queued, it automatically shuts 
down. This dynamic scaling optimizes cloud usage and helps reduce computational costs.

Since the EESD OpenSees distribution does not allow parallel execution, 
This Elastic Machine Group is composed of `c2d-highcpu-2` instances, each equipped with 2 virtual CPUs, since our goal is to 

```python
# Allocate an Elastic Machine Group on Google Cloud Platform
cloud_machine = inductiva.resources.ElasticMachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-2",
    spot=True,
    min_machines=1,
    max_machines=50)
```

> Learn more about the `ElasticMachineGroup` class [here](https://inductiva.ai/guides/how-it-works/machines/computational_resources/elasticgroup_class).

The elastic machine group 

e're using a cloud machine (`c2d-highcpu-4`) equipped with 4 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).


 In this example, we allocate an Elastic Machine Group composed of c2d-highcpu-2 instances. The group scales automatically from 1 to 50 machines depending on the workload. We chose this machine type, which has only 2 vCPUs, because this particular implementation of OpenSees is not parallelized. Using larger machines would not improve performance but would significantly increase costs.

### Section 2: Project Setup and Simulator Configuration
First, we specify the directory for the designated input files folder.

Next, we specify the OpenSees version to use. In this case, version 3.0.2, which includes the 3D macroelement, OrthotropicMembraneSection, and NoTensionSection3d commands, and supports the use of non-symmetric tangent stiffness in ND-zeroLength elements.

To ensure all simulation data is easily managed and to simplify post-processing and result downloads, it's best to associate all data with a project using the Inductiva Project class.

> Learn more about the Inductiva Project class [here](https://inductiva.ai/guides/scale-up/projects/projects).

```python
# Simulation input files
input_dir = r"C:/Path/To/Input/Files"

# Initialize the simulator
opensees = inductiva.simulators.OpenSees(
    interface="eesd",
    version="3.0.2")

# Organize all data into a Project
project_name = "B2_3P_Tutorial"
```

### Code Section 3: Computing Parameters for the Analysis
Next, we define the analysis range for the 30 records used in the IDA, along with the corresponding EQfactor values, which denote the level of intensity to scale the ground motion records.

We also define numerical damping, required for non-linear dynamic analyses. Rayleigh damping is specified in OpenSees via the Rayleigh command, which needs two parameters, alpha and beta, based on two fundamental frequencies of the structure. In this example, the Rayleigh damping is calibrated based on the 1<sup>st</sup> and 6<sup>th</sup> frequencies of the building. A damping ratio of 5% is assumed, typical for this type of analysis.

```python
# Variables of the analysis
# 30 bidirectional acceleration time-series to perform Incremental Dynamic Analysis
analysis_range = range(1,31)
# Scaling factor for the earthquake time-series
EQfactor_values = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

# Computing the parameters for a 5% fraction damping based on the 1st and 6th frequencies
damping_percentage = 0.05
Ti = 0.09970895596567399 # 1st vibration period
Tj = 0.04970502055069268 # 6th vibration period

alpha,beta = dp.damping(Ti, Tj, damping_percentage)
```

### Code Section 4: Initializing the Loop and `TemplateManager`
This section starts the loops over `analysis_range` and `EQfactor_values`. Each record is scaled to the specified intensity level, meaning the input files are rewritten at each iteration with updated parameters.

To automate this process, the `TemplateManager` tool is used. It overwrites the input files (with a .jinja extension) by replacing placeholders ({{ }}) with the appropriate parameter values at each step.

> Learn more about the `TemplateManager` [here](https://website-staging.inductiva.ai/guides/scale-up/parallel-simulations/templating).

```python
# Folder paths
# Including the .jinja files
input_files_template = os.path.join(input_dir, "inputFiles_template")
output_folder = os.path.join(input_dir, "outputFiles")
# Ensure output folder exists
os.makedirs(output_folder, exist_ok=True)

records_duration_file = os.path.join(input_dir, "records_duration.txt")
records_duration = np.loadtxt(records_duration_file, delimiter=' ')

# Loop over all combinations of ii and EQfactor
for ii in analysis_range:
    for EQfactor in EQfactor_values:

        print(f"Processing ii={ii}, EQfactor={EQfactor:.2f}")
        max_time = records_duration[ii-1]

        input_dir_folder = "Pedir ao daniel o script para eu ter este caminho completo, no print nao se ve"

        inductiva.TemplateManager.render_dir(
            source_dir=input_files_template,
            target_dir=input_dir_folder,
            overwrite=True,
            EQfactor=EQfactor,
            alpha=alpha,
            beta=beta,
            ii=ii,
            max_time=max_time)
```

### Code Section 5: Running Simulations and Defining Metadata
Within the loops, each simulation is executed by calling the batch file. Metadata is also defined for each run, capturing information essential for future post-processing.

```python
        batch_file_path = os.path.join(input_dir_folder, "Prototype_b2_3p_batch.tcl")
        shutil.copy(batch_file_path, input_dir)

        # Extract file name and extensision
        batch_file_name = os.path.basename(batch_file_path)
        name_batch, ext_batch = os.path.splitext(batch_file_name)
        
        # Run simulation
        task = opensees.run(
            input_dir=input_dir,
            sim_config_filename=batch_file_name,
            on=cloud_machine,
            project=project_name,
            resubmit_on_preemption=True)
        
        task.set_metadata({
            "factor": str(EQfactor),
            "Current_Analysis": str(ii),
            "Record_duration": str(max_time),
            "Project_name": project_name,
        })
```

### Code Section 6: Monitoring Simulations and Downloading Results
In the final stage, we wait for all simulations to complete. Once finished, the output data is automatically downloaded for further analysis, and the cloud machine group is terminated.

```python
inductiva.projects.Project(project_name).wait()
inductiva.projects.Project(project_name).download_outputs()
cloud_machine.terminate()
```

## Post-Processing and Analysis of Results
Once the simulation data is available, post-processing can be performed using tools of your choice—such as custom Python scripts.

In this tutorial, results are analyzed in terms of:
1. Capacity and fragility functions
2. Convergence of statistical moments (specifically, the median and dispersion)

Capacity curves are presented based on the 50th, 16th, and 84th percentiles of the IDA response, using the maximum global drift ratio as the engineering demand parameter (EDP). The 50th fractile represents the median capacity of the structure.

Next, the damage states of Immediate Occupancy (IO), Collapse Prevention (CP), and Near Collapse (NC) are identified. The IO limit is defined approximately at the end of the elastic branch. CP is identified when the slope of the capacity curve decreases to 20% of the initial elastic slope. Finally, NC is indicated by the appearance of the flatline.

<image2 and caption>

Fragility functions denote the likelihood of exceedance of a specified EDP value corresponding to each limit state, as shown in the figure below.

<image3 and caption>

The figures below illustrate the convergence of the median and dispersion of the fragility functions, based on the 30 observations for each limit state.

<image4/5 and caption>















# Run a Single Simulation
First, we will run a single OpenFAST simulation using the `5MW_OC4Semi_WSt_WavesWN` case. This should be straightforward as all 
the necessary input files are already prepared.

## Code Overview
The code required to run an OpenFAST simulation using the Inductiva API is always the same, we just need to adapt it for this use case, as shown below.

```python
import inductiva

# Allocate Cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
   provider="GCP",
   machine_type="n2-highcpu-2",
   spot=True)

# Initialize OpenFAST simulator
# Versions available: 3.5.2 and 4.0.2
openfast = inductiva.simulators.OpenFAST(
   version="4.0.2")

# Run simulation
task = openfast.run(
   input_dir="input_files",
   commands=[
       "openfast 5MW_OC4Semi_WSt_WavesWN/"
       "5MW_OC4Semi_WSt_WavesWN.fst"],
   on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as follows.

```
inductiva tasks info hk8fmb4dif05c5vf5saurocxp


Task status: Success


Timeline:
   Waiting for Input         at 20/02, 14:24:28      6.487 s
   In Queue                  at 20/02, 14:24:35      19.415 s
   Preparing to Compute      at 20/02, 14:24:54      1.842 s
   In Progress               at 20/02, 14:24:56      34.265 s
       └> 34.123 s        openfast 5MW_OC4Semi_WSt_WavesWN/5MW_OC4Semi_WSt_WavesWN.fst
   Finalizing                at 20/02, 14:25:30      1.085 s
   Success                   at 20/02, 14:25:31     


Data:
   Size of zipped output:    16.20 MB
   Size of unzipped output:  32.71 MB
   Number of output files:   83


Estimated computation cost (US$): 0.00013 US$
```

## Performance and Cost Analysis
Given that OpenFAST does not benefit from multiple CPU cores, we chose the `n2-highcpu-2` virtual machine (VM) with 2 virtual CPUs (equivalent to 1 physical core). 
This is one of the cheapest options on Google Cloud, costing just US$0.0081 per hour in spot mode.

To demonstrate that OpenFAST does not scale with the number of cores, we also ran the same simulation on a number of better machines. For a detailed breakdown, check out our [Benchmarks](../../benchmarks) section.

Here are the results:
| Machine       | vCPUs | Execution time | Estimated Cost (USD)|
|---------------|-----------------|----------------|------|
| n2-highcpu-2  | 2               |34.1s           |0.00013|
| n2-highcpu-4  | 4               |35.0s           |0.00031|
| n2-highcpu-8  | 8               |30.4s           |0.00034|
| n2-highcpu-16 | 16              |30.9s           |0.00065|
| n2-highcpu-32 | 32              |30.2s           |0.0012 | 

The execution time remains almost the same on all machines, regardless of the number of virtual CPUs.

Additionally, running the simulation on the `n2-highcpu-2` VM proves to be extremely cost-efficient, with a total cost of only 0.00013 US$.

In the next part of this tutorial, we'll take things to the next level by running dozens of OpenFAST simulations in parallel on Inductiva, demonstrating the true power of cloud-based scalability. Stay tuned!

```{banner_small}
:origin: openfast
```