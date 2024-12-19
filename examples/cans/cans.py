"""CaNS example"""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cans-input-example.zip",
    unzip=True)

# Initialize CaNS Simulator
# Check available versions with the cli command "inductiva simulators list"
cans = inductiva.simulators.CaNS(version="2.3.4")

# Run simulation with config files in the input directory
task = cans.run(input_dir=input_dir,
                sim_config_filename="input.nml",
                on=machine_group,
                n_vcpus=4)

task.wait()
task.download_outputs()

machine_group.terminate()
