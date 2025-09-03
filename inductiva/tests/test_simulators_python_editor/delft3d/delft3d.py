"""Delft3D example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import Delft3D
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "delft3d-input-example.zip",
    unzip=True)

# Initialize the Simulator
delft3d = Delft3D(version="6.04.00")

# List of commands to run
commands = ["mpirun -np 4 d_hydro.exe config_d_hydro.xml"]

task = delft3d.run( \
    input_dir=input_dir,
    commands=commands,
    on=machine)

task.wait(silent_mode=True)
machine.terminate()

print("=== Amazing! Your simulation has finished! ===")
