# Test Your Inductiva Setup
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick SFINCS simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code
To get started, copy the code below and paste it into a Python script.

When you run the script, all the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The SFINCS simulation will then be sent to a cloud machine for execution.

```python
"""SFINCS example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Download example configuration files from Inductiva storage
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "sfincs-input-example.zip",
    unzip=True)

# Initialize the Simulator
sfincs = inductiva.simulators.SFINCS( \
    version="2.2.1")

# Run simulation
task = sfincs.run( \
    input_dir=input_dir,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

<a href="https://console-dev.inductiva.ai/editor?simulator_name=sfincs" class="try-playground-button" target="_blank">
  <svg class="icon" xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
    <path d="M8 5v14l11-7z"/>
  </svg>
  Try it on our Python Editor, on any device
</a>

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal. If the task status shows **Success**, congratulations! You've successfully run an SFINCS simulation.

You're ready to start running simulations seamlessly!

```{banner_small}
:origin: sfincs
```

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
