"""WRF Simulation."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import WRF
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup( \
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/"
    "inductiva-api-demo-files/"
    "wrf-input-example.zip",
    unzip=True)

# Initialize the Simulator
wrf = WRF( \
    version="4.6.1")

# Run simulation
task = wrf.run( \
    input_dir=input_dir,
    init_commands=["./ideal.exe"],
    case_name="em_fire",
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
