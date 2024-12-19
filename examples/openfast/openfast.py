"""OpenFAST example."""
import inductiva
from datetime import timedelta

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    "c2-standard-4", spot=True, max_idle_time=timedelta(minutes=1))
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfast-input-example.zip",
    unzip=True)

# List of commands to run
commands = ["openfast IEA-15-240-RWT-Monopile.fst"]

# Initialize OpenFAST simulator
# Check available versions with the cli command "inductiva simulators list"
openfast = inductiva.simulators.OpenFAST(version="3.5.2")

# Run simulation
task = openfast.run(input_dir=input_dir, commands=commands, on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
