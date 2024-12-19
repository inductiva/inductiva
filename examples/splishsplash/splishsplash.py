"""SPlisHSPlasH example."""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip",
    unzip=True)

# Initialize SPlisHSPlasH Simulator
# Check available versions with the cli command "inductiva simulators list"
splishsplash = inductiva.simulators.SplishSplash(version="2.13.0")

task = splishsplash.run(input_dir=input_dir,
                        sim_config_filename="config.json",
                        on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
