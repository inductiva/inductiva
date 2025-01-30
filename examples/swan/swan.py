"""SWAN example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
swan = inductiva.simulators.SWAN()

# Run simulation
task = swan.run(input_dir="/path/to/my/swan/files",
                sim_config_filename="a11refr.swn",
                on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
