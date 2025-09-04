"""SNLSWAN example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import SNLSWAN
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "snlswan-input-example.zip", True)

# Initialize the Simulator
snlswan = SNLSWAN(version="2.2")

# Run simulation
task = snlswan.run( \
    input_dir=input_dir,
    sim_config_filename="input.swn",
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
