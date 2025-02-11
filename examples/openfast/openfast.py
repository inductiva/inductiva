"""OpenFAST example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    machine_type="c4-standard-4",
    provider="GCP")

# Initialize the Simulator
openfast = inductiva.simulators.OpenFAST()

my_openfast_command = [
    # List the OpenFAST commands you wish to execute
]

# Run simulation
task = openfast.run(input_dir="/path/to/my/openfast/files",
                    commands=my_openfast_command,
                    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
