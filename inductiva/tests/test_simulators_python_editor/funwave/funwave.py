"""FUNWAVE example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import FUNWAVE
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "funwave-input-example.zip",
    unzip=True)

# Initialize the Simulator
funwave = FUNWAVE( \
    version="3.6")

# Run simulation
task = funwave.run( \
    input_dir=input_dir,
    sim_config_filename="input.txt",
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
