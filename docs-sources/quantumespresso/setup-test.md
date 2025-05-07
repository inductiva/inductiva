# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick Quantum ESPRESSO simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code
To get started, copy the code below and paste it into a Python script.

When you run the script, all the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The Quantum ESPRESSO simulation will then be sent to a cloud machine for execution.

```python
"""Quantum ESPRESSO example."""
import inductiva
from inductiva.commands import MPIConfig, Command

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

mpi_config = MPIConfig( \
    version="4.1.6",
    np=2,
    use_hwthread_cpus=False)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "qe-input-example.zip",
    unzip=True)

# List of commands to run
commands = [
    Command("pw.x -i Al_local_pseudo.in", mpi_config=mpi_config),
    # openMP command should not be used with MPI
    "pw_openmp.x -i Al_qe_pseudo.in"
]

# Initialize the Simulator
qe = inductiva.simulators.QuantumEspresso( \
    version="7.4.1")

# Run simulation
task = qe.run( \
    input_dir=input_dir,
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal. If the task status shows **Success**, congratulations! You've successfully run a Quantum ESPRESSO simulation.

You're ready to start running simulations seamlessly!

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
