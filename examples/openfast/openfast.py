"""OpenFAST example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfast-input-example.zip",
    unzip=True)

# List of commands to run
commands = ["openfast IEA-15-240-RWT-Monopile.fst"]

# Initialize OpenFAST simulator
openfast = inductiva.simulators.OpenFAST()

# Run simulation
task = openfast.run(input_dir=input_dir, commands=commands, on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()

task.print_summary()
