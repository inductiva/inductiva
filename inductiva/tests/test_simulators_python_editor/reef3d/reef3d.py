"""Reef3D example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import REEF3D
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip",
    unzip=True)

# Initialize the Simulator
reef3d = REEF3D(version="24.02")

# Run simulation
task = reef3d.run( \
    input_dir=input_dir,
    on=machine)

task.wait(silent_mode=True)
machine.terminate()

print("=== Amazing! Your simulation has finished! ===")
