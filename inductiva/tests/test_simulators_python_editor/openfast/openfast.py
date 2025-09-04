"""OpenFAST example."""
from inductiva.resources.machine_groups import MachineGroup
from inductiva.simulators import OpenFAST
from inductiva.utils import download_from_url

# Instantiate machine group
machine = MachineGroup(machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfastv4.0.2-input-example.zip",
    unzip=True)

# List of commands to run
commands = ["openfast 5MW_OC4Semi_WSt_WavesWN/"
            "5MW_OC4Semi_WSt_WavesWN.fst"]

# Initialize the simulator
openfast = OpenFAST(version="4.0.2")

# Run simulation
task = openfast.run( \
    input_dir=input_dir,
    commands=commands,
    on=machine)

task.wait()
machine.terminate()

print("\n === Amazing! Your simulation has finished! ===")
