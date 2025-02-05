"""OpenFAST example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
openfast = inductiva.simulators.OpenFAST()

my_openfast_command = [
    # List the OpenFAST commands you wish to execute
]

# Run simulation
task = openfast.run(input_dir="/path/to/my/openfast/files",
                    commands=my_openfast_command,
                    on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()
