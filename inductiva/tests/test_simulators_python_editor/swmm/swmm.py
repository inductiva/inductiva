"""SWMM example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import SWMM
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-2")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "swmm-input-example.zip",
    unzip=True)

# Initialize the Simulator
swmm = SWMM(\
    version="5.2.4")

# List of commands to run
commands = [ \
    "runswmm model.inp model.rpt"]

# Run simulation
task = swmm.run(\
    input_dir=input_dir,
    commands=commands,
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
