"""XBeach example."""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))

machine_group.start()

# Download example configuration files from Inductiva storage
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "xbeach-input-example.zip",
    unzip=True)

# Initialize XBeach Simulator
# Check available versions with the cli command "inductiva simulators list"
xbeach = inductiva.simulators.XBeach(version="1.24")

# Run simulation with configuration files in the input directory
task = xbeach.run(input_dir=input_dir,
                  sim_config_filename="params.txt",
                  on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
