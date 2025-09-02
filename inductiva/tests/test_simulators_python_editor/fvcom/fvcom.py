"""FVCOM example"""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import FVCOM
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fvcom-input-example.zip",
    unzip=True)

# Initialize the Simulator
fvcom = FVCOM()

# Run simulation with config files in the input directory
task = fvcom.run( \
    input_dir=input_dir,
    working_dir="run/",
    case_name="tst",
    n_vcpus=1,
    debug=7,
    on=machine)

task.wait(silent_mode=True)
machine.terminate()

print("=== Amazing! Your simulation has finished! ===")
