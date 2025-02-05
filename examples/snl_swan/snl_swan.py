""" SNL SWAN example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
snl_swan = inductiva.simulators.SnlSwan()

# Run simulation with config files in the input directory
task = snl_swan.run(input_dir="/Path/to/My/Snl-Swan/Files",
                    sim_config_filename="my_config_file.swn",
                    on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()