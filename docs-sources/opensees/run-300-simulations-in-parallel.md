# IDAs on Inductiva: Run 300 Simulations in Parallel

*This tutorial was written by* [Daniel Caicedo](mailto:caicedod93@gmail.com) *in collaboration with the* **Inductiva Team**

<br>

**Incremental Dynamic Analysis (IDA)** is a computational method widely used in Performance-
Based Earthquake Engineering (PBEE) to evaluate the seismic performance of structures. This 
method generates response curves across a range of earthquake intensity levels, enabling a
probabilistic assessment of structural behaviour and associated seismic risk.

To date, the applicability of this technique for the case of unreinforced masonry constructions has been limited due to the **complexity of numerical models** and the **large computational burden** associated.

The first challenge is addressed by developing efficient structural models in OpenSees based on three-
dimensional macroelements to capture the in-plane (IP) and out-of-plane (OOP) mechanisms,
while also incorporating the effects of non-linear floor-to-wall connections and wall-to-wall
interlocking. 

The second, and more significant, challenge is tackled by leveraging the Inductiva
platform, which enables the execution of large-scale simulations through high-performance
computing (HPC) resources.

This tutorial showcases an example involving a two-storey unreinforced masonry building,
representative of the pre-code masonry building stock in metropolitan area of Lisbon. All the data
required for the development of numerical models was gathered in the scope of the [STAND4HERITAGE](https://stand4heritage.org) project (New STANDards for seismic assessment of built cultural HERITAGE).

<image2 and caption - ask Caicedo>

## Prerequisites
Download the required files [here](<link to input files>) and place them in a folder called `SimulationFiles`. 

-> Are these really prerequisites?
- Output files folder to save the results of the analysis.
- Records folder, where the 30 acceleration time-histories required for the IDAs are
located.
- Record duration file to be read during the execution of the script.

## Run 50 Simulations in Parallel
This section demonstrates how to run 50 Incremental Dynamic Analyses (IDAs) in parallel using an elastic group of cloud machines. The process is divided into six main code sections.

### Section 1: Allocating the Cloud Machine Group
We begin by allocating an elastic machine group to run the simulations. This group is configured with 
a minimum and maximum number of active machines. Machines are dynamically turned on or off
as demanded, ensuring efficient resource utilization.

Since some simulations may take longer than others, this elastic setup prevents idle waiting. When 
a machine finishes its assigned tasks and no further simulations are queued, it automatically shuts 
down. This dynamic scaling optimizes cloud usage and helps reduce computational costs.

```python
cloud_machine = inductiva.resources.ElasticMachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-2",
    spot=True,
    min_machines=1,
    max_machines=50)
```

>**Note**: In this example, we allocate an Elastic Machine Group composed of c2d-highcpu-2 instances. The group scales automatically from 1 to 50 machines depending on the workload. We chose this machine type, which has only 2 vCPUs, because this particular implementation of OpenSees is not parallelized. Using larger machines would not improve performance but would significantly increase costs.

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

<image2 and caption - ask Caicedo>

Fragility functions denote the likelihood of exceedance of a specified EDP value corresponding to each limit state, as shown in the figure below.

<image3 and caption - ask Caicedo>

The figures below illustrate the convergence of the median and dispersion of the fragility functions, based on the 30 observations for each limit state.

<image4/5 and caption - ask Caicedo>

## -> I think it's lacking a conclusion or an ending section