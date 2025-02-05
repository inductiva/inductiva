"""SWAN example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
swan = inductiva.simulators.SWAN(version="41.45")

# Run simulation with config files in the input directory
task = swan.run(input_dir="/path/to/my/swan/files",
                sim_config_filename="my_config_file.swn",
                on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()