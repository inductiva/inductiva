"""OpenFAST example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4")

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
task = openfast.run( \
    input_dir=input_dir,
    commands=commands,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
