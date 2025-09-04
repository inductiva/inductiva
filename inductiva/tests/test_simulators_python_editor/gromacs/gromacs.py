"""GROMACS example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import GROMACS
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "gromacs-input-example.zip",
    unzip=True)

commands = [
    "gmx solvate "
    "-cs tip4p -box 2.3 "
    "-o conf.gro -p topol.top",
]

# Initialize the Simulator
gromacs = GROMACS( \
    version="2022.2")

# Run simulation
task = gromacs.run( \
    input_dir=input_dir,
    commands=commands,
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
