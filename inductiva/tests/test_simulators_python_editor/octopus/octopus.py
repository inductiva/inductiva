"""Octopus example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import Octopus
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "octopus-input-example.zip",
    unzip=True)

# Initialize the Simulator
octopus = Octopus( \
    version="16.1")

# Run simulation
task = octopus.run( \
    input_dir=input_dir,
    commands=["octopus"],
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
