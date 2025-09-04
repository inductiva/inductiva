"""Quantum ESPRESSO example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import QuantumEspresso
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "qe-input-example.zip",
    unzip=True)

# List of commands to run
commands = ["pw_openmp.x -i Al_qe_pseudo.in"]

# Initialize simulator
qe = QuantumEspresso( \
    version="7.4.1")

# Run simulation
task = qe.run( \
    input_dir,
    commands=commands,
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
