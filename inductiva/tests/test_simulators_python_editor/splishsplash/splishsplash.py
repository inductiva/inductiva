"""SPlisHSPlasH example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import SplishSplash
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip",
    unzip=True)

# Initialize the Simulator
splishsplash = SplishSplash()

# Run simulation
task = splishsplash.run( \
    input_dir=input_dir,
    sim_config_filename="config.json",
    on=machine)

task.wait(silent_mode=True)
machine.terminate()

print("=== Amazing! Your simulation has finished! ===")