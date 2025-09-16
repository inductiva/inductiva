# Kickstart Your Simulation Workflow

## Prerequisites
To follow this tutorial, either as a hands-on walkthrough or as a reference for your own implementation, please begin by downloading 
the simulation files [here](https://storage.googleapis.com/inductiva-api-demo-files/opensees-tutorials/IDA-at-scale.zip).

The files are organized into three folders to support both pre- and post-processing:
- `inputFiles_template` – contains the FEM model definition and analysis scripts, including the dynamic analysis template 
(`.jinja` file) that will be modified within the script’s loops
- `outputFiles` – the directory where the analysis results will be saved
- `Records` – contains the 30 acceleration time-histories required for the IDAs

A record duration file, `records_duration.txt`, is also included and will be read during the execution of the simulation script.

## Code Overview
The Python code provided below is designed to run OpenSees simulations on Inductiva’s cloud infrastructure. It is fully adaptable 
and can be customized for virtually any OpenSees use case.

In this example, we demonstrate how to run **300 Incremental Dynamic Analyses (IDAs) in parallel**, showcasing a typical 
high-throughput simulation workflow using the **EESD OpenSees distribution**.

The code is divided into seven key sections:
1. Defining the Damping Function
2. Allocating the Cloud Machine Group
3. Setting Up the Project and Configuring the Simulator
4. Computing Analysis Parameters
5. Initializing the Simulation Loop and `TemplateManager`
6. Running the Simulations and Assigning Metadata
7. Monitoring Progress and Downloading Results

We’ll guide you through each of these sections to help you understand both the structure and flexibility of the workflow.

Copy the code and save it as a Python file - your simulation script - for execution after reviewing the following sections.

```python
import numpy as np
import os
import shutil
import inductiva
import math

# Defining the Damping Function
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

# Allocating the Cloud Machine Group
cloud_machine = inductiva.resources.ElasticMachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-2",
    spot=True,
    min_machines=1,
    max_machines=50)

# Setting up the project and configuring the simulator
input_dir = r"/Path/to/Tutorial/Files"

opensees = inductiva.simulators.OpenSees(
    interface="eesd",
    version="3.0.2")

project_name = "B2_3P_Tutorial"

# Computing Analysis Parameters
analysis_range = range(1,31) 

EQfactor_values = [
    0.05, 0.10, 0.15,
    0.20, 0.25, 0.30,
    0.35, 0.40, 0.45,
    0.50]

damping_percentage = 0.05
Ti = 0.09970895596567399
Tj = 0.04970502055069268

alpha,beta = damping(Ti, Tj, damping_percentage)

# Initializing the simulation loop and TemplateManager
input_files_template = os.path.join(
    input_dir, "inputFiles_template")
output_folder = os.path.join(
    input_dir, "outputFiles")

os.makedirs(output_folder, exist_ok=True)

records_duration_file = os.path.join(
    input_dir, "records_duration.txt")
records_duration = np.loadtxt(
    records_duration_file, delimiter=' ')

for ii in analysis_range:
    for EQfactor in EQfactor_values:

        print(f"Processing ii={ii}, EQfactor={EQfactor:.2f}")
        max_time = records_duration[ii-1]

        #Directory where the rendered simulation files will be placed
        input_dir_folder = r"/Path/to/Tutorial/Files/inputFiles"

        inductiva.TemplateManager.render_dir(
            source_dir=input_files_template,
            target_dir=input_dir_folder,
            overwrite=True,
            EQfactor=EQfactor,
            alpha=alpha,
            beta=beta,
            ii=ii,
            max_time=max_time)

        # Running the simulations and assigning metadata
        batch_file_path = os.path.join(
            input_dir_folder, "Prototype_b2_3p_batch.tcl")
        shutil.copy(
            batch_file_path, input_dir)

        batch_file_name = os.path.basename(
            batch_file_path)
        name_batch, ext_batch = os.path.splitext(batch_file_name)
                
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
            "Project_name": project_name})

# Monitoring progress and downloading results
inductiva.projects.Project(project_name).wait()
inductiva.projects.Project(project_name).download_outputs()
cloud_machine.terminate()
```
