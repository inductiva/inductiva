# Test Your Inductiva Setup
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick OpenSees simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code
To get started, copy the code below and paste it into a Python script.

When you run the script, all the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The OpenSees simulation will then be sent to a cloud machine for execution.

```python
"""OpenSees example."""
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openseespy-input-example.zip",
    unzip=True)

# Initialize the Simulator
opensees = inductiva.simulators.OpenSees( \
    interface="python",
    version="3.7.1")

# Run simulation
task = opensees.run( \
    input_dir=input_dir,
    sim_config_filename="example_mpi_paralleltruss_explicit.py",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal. If the task status shows **Success**, congratulations! You've successfully run an OpenSees simulation.

You're ready to start running simulations seamlessly!

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!

```{banner_small}
:origin: opensees
```
