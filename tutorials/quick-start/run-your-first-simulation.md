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

<div style="border: 2px solid #007BFF; padding: 15px; background-color: #F0F8FF; border-radius: 8px; margin: 20px 0;">
  <h3 style="margin-top: 0;">🎥 Prefer a visual guide?</h3>
  <p>Follow along with our <strong>video walkthrough</strong> led by our CEO, where key details, features, and step-by-step instructions are explained in depth. It's a quick and engaging way to get started in no time!</p>
</div>

Let’s go!

## Prerequisites

Before starting, ensure the following:

* **Python Installed:** Make sure Python (>=3.8) is installed on your system. Check our [System Prep guide](https://docs.inductiva.ai/en/latest/preinstallation/system/system-requirements.html#select-your-os) for more information.
* **Inductiva API Setup:** 

    - Install the Inductiva library using `pip install inductiva`.
    - Authenticate your API key by [following these instructions](https://docs.inductiva.ai/en/latest/preinstallation/system/system-requirements.html#installing-the-inductiva-python-client-and-authenticating).


## Step 1: Prepare Your Project Directory

### Get Your Simulation Artifacts

1.	**Locate the Input Files**
   
    Download the necessary input files from the [Galveston Island dataset](https://data.griidc.org/data/HI.x833.000:0001#individual-files).

2.	**Navigate to the Correct Folder**

    Find the files under: *Files >> XBeach_Model_Runs >> Beach_Nourish_Only >> Input.*

3. Download all the files into a folder named **Beach_Nourish_Only** on your local machine.

### Adjust Simulation Parameters

1. **Open `params.txt`and modify the parameters**

    Ensure compatibility with XBeach v10+ by adding the following line after the header:

    `single_dir = 0`

    Reduce simulation time by changing the tstop value:

    `tstop = 34560`

### Create the Python Script (run.py)

Open your preferred IDE (e.g., VS Code, PyCharm) and save a new file as `run.py` in your project directory.

<div style="border: 2px solid #28A745; padding: 15px; background-color: #E9FBEA; border-radius: 8px; margin: 20px 0;">
  <h3 style="margin-top: 0;">💡 Pro Tip: Quick Creation of `run.py`</h3>
  <p>Instead of creating the file manually, you can use a shortcut:  
     - On <strong>Mac/Linux</strong>, type:  
       <code>touch run.py</code>  
     - On <strong>Windows</strong> (PowerShell), type:  
       <code>New-Item -ItemType File -Name "run.py"</code>  
  </p>
  <p>This will instantly create the `run.py` file in your current directory.</p>
</div>

### Organise Your Files

Ensure your directory structure looks like this:

```
xbeach/
|-- Beach_Nourish_Only/
    |-- README.txt
    |-- bed.dep
    |-- bedfricfile.txt
    ...
|-- run.py
```

* The `beach_nourish_only` folder contains all simulation artifacts.

* The`run.py` file which will contain your script to execute the simulation.

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

### Execute the Script

Run the script from your terminal:

```bash
python3 run.py
```
![Live simulation progress](../_static/simulation_progress.gif "Live simulation progress displayed in the terminal, showing machine setup, file uploads, and real-time updates. You can monitor the logs here or follow along with our video walkthrough for a detailed breakdown.")

The script will:

- Authenticate your API key.
- Start the specified cloud machine in spot mode.
- Zip and upload simulation files to Inductiva servers
- Execute the simulation on the cloud.
- Wait for completion and terminate the machine.

### Monitor the Simulation

**Command Line:** Use the inductiva logs command to view simulation logs:

```bash
inductiva logs <TASK_ID>
```

**Web Console:** Log in to the Inductiva Console to track simulation progress and view outputs.

## Step 4: Accessing Outputs

- Download outputs via the Web Console or CLI.
- Partial outputs are saved if the simulation is interrupted.
- Use the outputs as inputs for subsequent simulations if needed.

## Good to Know! 

### Cost Management

- Spot mode significantly reduces costs (~85 cents for this example).
- Idle machines are automatically shut down after 30 minutes to prevent unexpected charges.

## Conclusion

By completing this tutorial, you’ve successfully run an advanced XBeach simulation using the Inductiva API. Explore other use cases and simulators available on the Inductiva Tutorials Page. If you encounter issues, refer to the documentation or contact support.

Happy simulating!