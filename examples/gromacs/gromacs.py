"""GROMACS example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

my_gmx_command = [
    # here you list the gromacs command you wish to execute
]

# Initialize the Simulator
gromacs = inductiva.simulators.GROMACS()

task = gromacs.run(input_dir="path/to/my/gromacs/files",
                   commands=my_gmx_command,
                   on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
