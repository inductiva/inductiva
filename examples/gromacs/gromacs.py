"""GROMACS example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-180")

my_gmx_command = [
    # List the GROMACS commands you wish to execute
]

# Initialize the Simulator
gromacs = inductiva.simulators.GROMACS()

# Run simulation with config files in the input directory
task = gromacs.run(input_dir="path/to/my/gromacs/files",
                   commands=my_gmx_command,
                   on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
