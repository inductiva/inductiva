"""GprMax example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import GprMax
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "gprmax-input-example.zip",
    unzip=True)

# Initialize the Simulator
gprmax = GprMax(\
    version="3.1.7")

# List of commands to run
commands = [ \
    "python -m gprMax antenna_like_MALA_1200_fs.in"]

# Run simulation
task = gprmax.run(\
    input_dir=input_dir,
    commands=commands,
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
