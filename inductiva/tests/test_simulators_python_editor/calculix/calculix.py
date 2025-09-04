"""Calculix example"""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import Calculix
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "calculix-input-example.zip",
    unzip=True)

# Initialize the Simulator
calculix = Calculix( \
    version="2.22")

# Run simulation
task = calculix.run( \
    input_dir=input_dir,
    sim_config_filename="hueeber3.inp",
    on=machine,
    n_vcpus=4)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
