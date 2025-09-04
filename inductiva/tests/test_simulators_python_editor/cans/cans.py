"""CaNS example"""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import CaNS
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cansv2.4.0-input-example.zip",
    unzip=True)

# Initialize the Simulator
cans = CaNS(version="2.4.0")

# Run simulation
task = cans.run( \
    input_dir=input_dir,
    sim_config_filename="input.nml",
    on=machine,
    n_vcpus=4)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
