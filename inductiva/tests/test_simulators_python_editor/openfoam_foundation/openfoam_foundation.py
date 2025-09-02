"""OpenFOAM Foundation example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import OpenFOAM
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-input-example.zip",
    unzip=True)

# Initialize the Simulator
openfoam = OpenFOAM(distribution="foundation", version="8")

# Run simulation
task = openfoam.run( \
    input_dir=input_dir,
    shell_script="./Allrun",
    on=machine)

task.wait(silent_mode=True)
machine.terminate()

print("=== Amazing! Your simulation has finished! ===")
