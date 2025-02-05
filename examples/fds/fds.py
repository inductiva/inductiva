"""FDS example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
fds = inductiva.simulators.FDS()

# Run simulation with config files in the input directory
task = fds.run(input_dir="path/to/my/fds/files",
               sim_config_filename="my_config_file.fds",
               on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()