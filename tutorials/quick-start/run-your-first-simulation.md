# Run Your First XBeach Simulation Example with Inductiva 

Welcome to Inductiva! Now that you’ve completed the [onboarding process](https://console.inductiva.ai/), 
it’s time to dive in and see how easy it is to run your first simulation!

This tutorial demonstrates how to run a simulation using the XBeach simulator 
through the Inductiva API. We will use a real-world dataset from the Galveston Island, Texas 
study on beach and dune enhancement scenarios [(dataset available here)](https://data.griidc.org/data/HI.x833.000:0001). 

By following this step-by-step guide, you will:

* Learn how to set up and execute a simulation using the Inductiva API.
* Use and select high-performance cloud machines.
* Access and manage simulation outputs via CLI and our [Web Console](https://console.inductiva.ai/).

🎥 **Prefer a visual guide?** Watch our video walkthrough to follow along step by step and get started in no time!

Let’s go!

## Prerequisites

Before starting, ensure the following:

* **Python Installed:** Make sure Python (>=3.8) is installed on your system.
* **Inductiva API Setup:** 

    - Install the Inductiva library using `pip install inductiva`.
    - Authenticate your API key by [following these instructions](https://docs.inductiva.ai/en/latest/preinstallation/system/system-requirements.html#installing-the-inductiva-python-client-and-authenticating).

* **Simulation Files:**
    - Locate the necessary input files from the [Galveston Island dataset](https://data.griidc.org/data/HI.x833.000:0001#individual-files).
    - Under Files >> XBeach_Model_Runs >> Beach_Nourish_Only >> Input, download all files into a a Beach_Nourish_Only folder.
    - Organize them in a directory structure as shown below.

## Step 1: Set Up Your Environment

### Folder Structure

Organize your files as follows:

```
xbeach/
|-- Beach_Nourish_Only/
    |-- README.txt
    |-- bed.dep
    |-- bedfricfile.txt
    ...
|-- run.py

```

* The `beach_nourish_only` folder contains all simulation artifacts downloaded from [GRIIDC] (https://data.griidc.org/data/HI.x833.000:0001#individual-files).
* An empty `run.py` file which will contain your script to execute the simulation. *You can create it by saving a new python file using your preferred IDE or simply type `touch run.py` on the command line if on Mac or Linux.*

## Step 2: Write the Python Script

In the empty `run.py` file you've created, add the following code, save, and close!

```py
# Start by importing the required Inductiva library
import inductiva

# Define the machine group
machine_group = inductiva.resources.MachineGroup(
    machine_type="c3d-highcpu-90",
    spot=True, # Enables cost-saving spot mode
    data_disk_gb=20) 
machine_group.start() # Start the machine group

# Specify the input directory
input_dir = "Beach_Nourish_Only"

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach()

# Run simulation with config files in the input directory
task = xbeach.run(
    input_dir=input_dir,
    sim_config_filename="params.txt",
    n_vcpus=90,
    on=machine_group)

# task.wait() is a blocking call and will only return when the simulation
# ends. However, you can close your terminal without interrupting the 
# simulation and use Inductiva CLI (Command Line Interface) tools to
# check the status of the simulation from another terminal.
task.wait()

# Clean up resources after the simulation
machine_group.terminate()

# Let's get a small summary of the run.
task.print_summary()

```

## Step 3: Run the Simulation


## Common Issues and Troubleshooting

Problem: Simulation not running?
Solution: Ensure all input parameters are defined correctly.

Problem: Output looks strange?
Solution: Check your initial conditions or domain size.

## Next Steps

Check out these resources:

- Inductiva Documentation
- Simulation Gallery