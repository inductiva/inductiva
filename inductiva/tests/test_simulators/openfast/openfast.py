"""OpenFAST example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfastv4.0.2-input-example.zip",
    unzip=True)

# List of commands to run
commands = ["openfast IEA-15-240-RWT-Monopile.fst"]

# Initialize OpenFAST simulator
openfast = inductiva.simulators.OpenFAST(
    version="4.0.2")

# Run simulation
task = openfast.run(input_dir=input_dir, commands=commands, on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
