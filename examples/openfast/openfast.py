"""OpenFAST example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c4-standard-4")

# Initialize the Simulator
openfast = inductiva.simulators.OpenFAST()

# Specify the OpenFAST commands you want to run, separated by commas
# Example using the 'openfast' command
openfast_commands = ["openfast my_file.fst"]

# Run simulation
task = openfast.run( \
    input_dir="/path/to/my/openfast/files",
    commands=openfast_commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
