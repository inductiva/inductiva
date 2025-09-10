"""Opentelemac example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import OpenTelemac
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "opentelemac-input-example.zip",
    unzip=True)

# Initialize the Simulator
opentelemac = OpenTelemac( \
    version="8p4r0")

# List of commands to run
commands = [ \
    "telemac2d.py t2d_malpasset-fine.cas --ncsize=16"]

# Run simulation
task = opentelemac.run( \
    input_dir=input_dir,
    commands=commands,
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
