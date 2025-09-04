"""SWAN example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import SWAN
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "swan-input-example.zip", True)

# Initialize the Simulator
swan = SWAN( \
    version="41.45")

# Run simulation
task = swan.run( \
    input_dir=input_dir,
    sim_config_filename="a11refr.swn",
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
