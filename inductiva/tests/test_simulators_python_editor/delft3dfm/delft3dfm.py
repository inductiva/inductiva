"""Delft3DFM example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import Delft3DFM
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "delft3dfm-input-example.zip",
    unzip=True)

# Initialize the Simulator
delft3dfm = Delft3DFM( \
    version="2024.03")

# List of commands to run
commands = [ \
    "dflowfm --autostartstop f34.mdu"]

task = delft3dfm.run( \
    input_dir=input_dir,
    commands=commands,
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
