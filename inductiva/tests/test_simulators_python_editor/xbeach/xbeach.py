"""XBeach example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import XBeach
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Download example configuration files from Inductiva storage
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "xbeach-input-example.zip",
    unzip=True)

# Initialize the Simulator
xbeach = XBeach(version="1.24")

# Run simulation
task = xbeach.run( \
    input_dir=input_dir,
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
