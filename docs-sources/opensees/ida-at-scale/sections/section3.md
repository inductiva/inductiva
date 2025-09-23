# Part 2: Running the Simulations on the Cloud
With the case fully defined, we now turn to executing the simulations at scale using the Inductiva cloud platform. This includes:
* Allocating compute resources
* Configuring the simulator
* Running the simulation loop
* Monitoring progress and retrieving results

## 1. Allocating the Cloud Machine Group
We allocate an **Elastic Machine Group** to run the simulations. This group is configured with a minimum and maximum 
number of machines, which automatically turn on or off based on workload, ensuring efficient use of resources.

Since some simulations may take longer than others, this elastic setup avoids idle waiting. Once a machine completes its assigned tasks and no further simulations are queued, it automatically shuts down. This dynamic scaling optimizes cloud usage and helps reduce computational costs.

Because the **EESD OpenSees distribution** runs on a **single thread**, it cannot take full advantage of multi-core CPUs. To maximize throughput, it’s more efficient to run multiple small machines, such as `c2d-highcpu-2` instances with 2 virtual CPUs, concurrently. This approach enables high parallelism while keeping resource usage and cost under control.

```python
# Allocating the Cloud Machine Group
cloud_machine = inductiva.resources.ElasticMachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-2",
    spot=True,
    min_machines=1,
    max_machines=50)
```

> Learn more about the `ElasticMachineGroup` class [here](https://inductiva.ai/guides/how-it-works/machines/computational_resources/elasticgroup_class).

## 2. Simulator Configuration
Once the cloud machines are allocated, we move on to defining the simulation inputs and selecting the simulator.

We start by specifying the directory that contains the input files for the analyses.

Next, we specify the EESD OpenSees version to use. In this case, version 3.0.2, which includes the 3D macroelement, `OrthotropicMembraneSection`, and `NoTensionSection3d` commands, as well as support for non-symmetric tangent stiffness in ND-zeroLength elements.

```python
# Initialize the Simulator
opensees = inductiva.simulators.OpenSees(
    interface="eesd",
    version="3.0.2")
```

## 3. Initializing the Simulation Loop and Running the Simulations
This section goes through every combination of analysis case (`analysis_range`) and scaling factor (`EQfactor_values`). Each record is scaled to the specified intensity level and then 
the input files are rewritten at each iteration with updated parameters.

To automate this process, the `TemplateManager` tool is used. It overwrites the `.jinja` input files by replacing placeholders ({{ }}) with the appropriate parameter values at each step.

Within the loops, each simulation is executed by calling the batch file. Metadata is also defined for each run, capturing information essential for future post-processing.

```python
for ii in analysis_range:
    for EQfactor in EQfactor_values:

        print(f"Processing ii={ii}, EQfactor={EQfactor:.2f}")
        max_time = records_duration[ii - 1]

        #Render the templates into simulation files with
        #EQfactor, alpha, beta, ii and max_time values replaced
        inductiva.TemplateManager.render_dir(source_dir=input_files_template,
                                             target_dir=simulation_files,
                                             overwrite=True,
                                             EQfactor=EQfactor,
                                             alpha=alpha,
                                             beta=beta,
                                             ii=ii,
                                             max_time=max_time)
        # Running the Simulations
        task = opensees.run(input_dir=simulation_files,
                            sim_config_filename=tcl_file,
                            on=cloud_machine,
                            project=project_name,
                            resubmit_on_preemption=True)

        #Assigning Metadata
        task.set_metadata({
            "factor": str(EQfactor),
            "Current_Analysis": str(ii),
            "Record_duration": str(max_time),
            "Project_name": project_name
        })
```

> Learn more about the `TemplateManager` class [here](https://inductiva.ai/guides/scale-up/parallel-simulations/templating).

## 4. Monitoring Progress and Downloading Results
In the final stage, we wait for all simulations to complete. Once finished, the output data is automatically downloaded for 
further analysis, and the cloud machine group is terminated.

```python
# Monitoring Progress and Downloading Results
inductiva.projects.Project(project_name).wait()
inductiva.projects.Project(project_name).download_outputs()
cloud_machine.terminate()
```

Now let’s move on to the final tutorial section to analyze the results.