"""FDS example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import FDS
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "fds-input-example.zip",
    unzip=True)

# Initialize the Simulator
fds = FDS( \
    version="6.10.1")

# Run simulation
task = fds.run( \
    input_dir=input_dir,
    sim_config_filename="mccaffrey.fds",
    n_vcpus=1,
    on=machine)

task.wait()
machine.terminate()
task.print_summary()

print("\n === Amazing! Your simulation has finished! ===")
