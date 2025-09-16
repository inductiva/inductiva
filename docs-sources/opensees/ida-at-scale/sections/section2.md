# Run 300 Simulations in Parallel
Let's walk through each of the seven sections of the simulation script.

## Section 1: Defining the Damping Function
We begin by defining a **damping function** that calculates the coefficients α and β, which are essential for modeling energy dissipation 
in structural systems.

This function takes as input the periods of two vibration modes and a target damping ratio, then returns coefficients that ensure 
consistent damping across a range of frequencies.

These coefficients are fundamental in OpenSees simulations, as they enable a realistic representation of how structures lose 
vibrational energy under dynamic loading.

```python
def damping(Ti, Tj, ksi):
    """
    Calculate damping coefficients alpha and beta.

    Parameters:
        Ti (float): Period of the first mode.
        Tj (float): Period of the second mode.
        ksi (float): Damping ratio.

    Returns:
        tuple: (alpha, beta) damping coefficients.
    """
    fi = 1 / Ti
    fj = 1 / Tj

    wi = 2 * math.pi * fi
    wj = 2 * math.pi * fj

    alpha = ksi * 2 * wi * wj / (wi + wj)
    beta = ksi * 2 / (wi + wj)

    return alpha, beta
```

## Section 2: Allocating the Cloud Machine Group
Next, we allocate an **Elastic Machine Group** to run the simulations. This group is configured with a minimum and maximum 
number of machines, which automatically turn on or off based on workload, ensuring efficient use of resources.

Since some simulations may take longer than others, this elastic setup avoids idle waiting. Once a machine completes its assigned 
tasks and no further simulations are queued, it automatically shuts down. This dynamic scaling optimizes cloud usage and helps reduce computational costs.

Because the **EESD OpenSees distribution** runs on a **single thread**, it cannot take full advantage of multi-core CPUs. 
To maximize throughput, it’s more efficient to run multiple small machines, such as `c2d-highcpu-2` instances with 2 virtual 
CPUs, concurrently. This approach enables high parallelism while keeping resource usage and cost under control.

```python
# Allocate an Elastic Machine Group on Google Cloud Platform
cloud_machine = inductiva.resources.ElasticMachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-2",
    spot=True,
    min_machines=1,
    max_machines=300)
```

> Learn more about the `ElasticMachineGroup` class [here](https://inductiva.ai/guides/how-it-works/machines/computational_resources/elasticgroup_class).

## Section 3: Project Setup and Simulator Configuration
Once the cloud machines are allocated, we move on to defining the simulation inputs and selecting the simulator.

We start by specifying the directory that contains the input files for the analyses.

Next, we specify the EESD OpenSees version to use. In this case, version 3.0.2, which includes the 3D macroelement, `OrthotropicMembraneSection`, and `NoTensionSection3d` commands, as well as support for non-symmetric tangent stiffness in 
ND-zeroLength elements.

To keep simulation data organized and make result handling easier, we group everything into a project using the Project class.

```python
# Input files path
input_dir = r"/Path/to/Tutorial/Files"

# Initialize the Simulator
opensees = inductiva.simulators.OpenSees(
    interface="eesd",
    version="3.0.2")

# Organize all data into a Project
project_name = "B2_3P_Tutorial"
```

> Learn more about the `Project` class [here](https://inductiva.ai/guides/scale-up/projects/projects).

## Section 4: Computing Analysis Parameters
Next, we define the analysis range for the 30 ground motion records used in the IDA, along with their corresponding EQfactor values, which represent the intensity scaling levels for each record. For every one of the 30 earthquake records, sourced from events around the world, we apply 10 different intensity levels.

We also define numerical damping, required for non-linear dynamic analyses. Rayleigh damping is specified in OpenSees via the 
Rayleigh command, which needs two parameters, alpha and beta, based on two fundamental frequencies of the structure. In this example, 
the Rayleigh damping is calibrated based on the 1<sup>st</sup> and 6<sup>th</sup> frequencies of the building. A damping ratio of 5% 
is assumed, typical for this type of analysis.

```python
# Analysis variables
# 30 bidirectional acceleration time-series to perform IDA
analysis_range = range(1,31)

# Scaling factor for the earthquake time-series
EQfactor_values = [
    0.05, 0.10, 0.15,
    0.20, 0.25, 0.30,
    0.35, 0.40, 0.45,
    0.50]

# Computing the parameters for a 5% fraction
# damping based on the 1st and 6th frequencies
damping_percentage = 0.05
Ti = 0.09970895596567399 # 1st vibration period
Tj = 0.04970502055069268 # 6th vibration period

alpha,beta = dp.damping(Ti, Tj, damping_percentage)
```

## Section 5: Initializing the Simulation Loop and `TemplateManager`
This section starts the loops over `analysis_range` and `EQfactor_values`. Each record is scaled to the specified intensity level, 
meaning the input files are rewritten at each iteration with updated parameters.

To automate this process, the `TemplateManager` tool is used. It overwrites the `.jinja` input files by replacing placeholders ({{ }}) 
with the appropriate parameter values at each step.

```python
# Including the .jinja files
input_files_template = os.path.join(
    input_dir, "inputFiles_template")
output_folder = os.path.join(
    input_dir, "outputFiles")

# Ensure output folder exists
os.makedirs(
    output_folder, exist_ok=True)

records_duration_file = os.path.join(
    input_dir, "records_duration.txt")
records_duration = np.loadtxt(
    records_duration_file, delimiter=' ')

# Initializing the Simulation Loop
# Loop over all combinations of ii and EQfactor
for ii in analysis_range:
    for EQfactor in EQfactor_values:

        print(f"Processing ii={ii}, EQfactor={EQfactor:.2f}")
        max_time = records_duration[ii-1]

        #Place where the rendered simulation files will be placed
        input_dir_folder = r"/Path/to/Tutorial/Files/inputFiles"

        #Render the templates into simulation files with
        #EQfactor, alpha, beta and ii values replaced
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

> Learn more about the `TemplateManager` class [here](https://inductiva.ai/guides/scale-up/parallel-simulations/templating).

## Section 6: Running the Simulations and Assigning Metadata
Within the loops, each simulation is executed by calling the batch file. Metadata is also defined for each run, capturing 
information essential for future post-processing.

```python
        batch_file_path = os.path.join(
            input_dir_folder, "Prototype_b2_3p_batch.tcl")
        shutil.copy(
            batch_file_path, input_dir)

        # Extract file name and extension
        batch_file_name = os.path.basename(
            batch_file_path)
        name_batch, ext_batch = os.path.splitext(
            batch_file_name)
        
        # Running the Simulations
        task = opensees.run(
            input_dir=input_dir,
            sim_config_filename=batch_file_name,
            on=cloud_machine,
            project=project_name,
            resubmit_on_preemption=True)
        
        #Assigning Metadata
        task.set_metadata({
            "factor": str(EQfactor),
            "Current_Analysis": str(ii),
            "Record_duration": str(max_time),
            "Project_name": project_name,
        })
```

## Section 7: Monitoring Progress and Downloading Results
In the final stage, we wait for all simulations to complete. Once finished, the output data is automatically downloaded for 
further analysis, and the cloud machine group is terminated.

```python
# Monitoring Progress and Downloading Results
inductiva.projects.Project(project_name).wait()
inductiva.projects.Project(project_name).download_outputs()
cloud_machine.terminate()
```

Now that we’ve covered each of the seven code sections in detail, let’s move on to the final tutorial section to analyze the results.