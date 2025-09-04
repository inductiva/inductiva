"""CM1 example"""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import CM1
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "cm1-input-example.zip",
    unzip=True)

# Initialize the Simulator
cm1 = CM1( \
    version="21.1")

# Run simulation
task = cm1.run( \
    input_dir=input_dir,
    sim_config_filename="namelist.input",
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
