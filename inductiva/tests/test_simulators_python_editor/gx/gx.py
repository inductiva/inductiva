"""Gx example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import GX
from inductiva.utils import download_from_url

# Instantiate machine group
gpu_machine = MachineGroup( \
    machine_type="g2-standard-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "gx-input-example.zip",
    unzip=True)

# Initialize the Simulator
gx = GX( \
    version="11-2024")

# Run simulation
task = gx.run( \
    input_dir=input_dir,
    sim_config_filename="itg_w7x_adiabatic_electrons.in",
    on=gpu_machine)

task.wait()
gpu_machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
