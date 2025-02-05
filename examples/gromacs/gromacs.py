"""GROMACS example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

my_gmx_command = [
    # List the GROMACS commands you wish to execute
]

# Initialize the Simulator
gromacs = inductiva.simulators.GROMACS()

# Run simulation
task = gromacs.run(input_dir="path/to/my/gromacs/files",
                   commands=my_gmx_command,
                   on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()