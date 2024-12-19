"""FDS example."""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fds-input-example.zip",
    unzip=True)

# Initialize FDS Simulator
# Check available versions with the cli command "inductiva simulators list"
fds = inductiva.simulators.FDS(version="6.8")

task = fds.run(input_dir=input_dir,
               sim_config_filename="mccaffrey.fds",
               n_vcpus=1,
               on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
