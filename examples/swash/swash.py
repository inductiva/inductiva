"""SWASH example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
swash = inductiva.simulators.SWASH(version="10.05")

# Run simulation with config files in the input directory
task = swash.run(input_dir="/path/to/my/swash/files",
                 sim_config_filename="input.sws",
                 on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
