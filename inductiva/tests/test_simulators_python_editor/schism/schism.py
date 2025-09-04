"""SCHISM example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import SCHISM
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "schism-input-example.zip", True)

# Initialize the Simulator
schism = SCHISM( \
    version="5.11.0")

# Run simulation
task = schism.run( \
    input_dir=input_dir,
    n_vcpus=3,
    num_scribes=2,
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
