"""SWASH example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
swash = inductiva.simulators.SWASH(version="10.05")

# Run simulation with config files in the input directory
task = swash.run(input_dir="/path/to/my/swash/files",
                 sim_config_filename="my_config_file.sws",
                 on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()