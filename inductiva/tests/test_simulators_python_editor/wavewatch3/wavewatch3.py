"""WAVEWATCH III example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import WaveWatch3
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4", spot=True)

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "wavewatch3-input-example.zip",
    unzip=True)

# Initialize the Simulator
wavewatch3 = WaveWatch3()

# Define the commands to run the simulation
commands = [ \
    "ww3_grid",
    "ww3_prep",
    "ww3_shel"]

# Run the simulation
task = wavewatch3.run( \
    input_dir=input_dir,
    custom_switch="switch_PR3_UQ_MPI",
    commands=commands,
    on=machine)

task.wait(silent_mode=True)
machine.terminate()

print("=== Amazing! Your simulation has finished! ===")
