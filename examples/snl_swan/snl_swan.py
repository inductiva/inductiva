""" SNL SWAN example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")


# Initialize the Simulator
snl_swan = inductiva.simulators.SnlSwan()

# Run simulation with config files in the input directory
# Uses swanrun by default
task = snl_swan.run(input_dir="/Path/to/My/Snl-Swan/Files",
                sim_config_filename="INPUT.swn",
                on=machine_group)

# Uses swanrun
# task = snl_swan.run(input_dir="/Path/to/My/Snl-Swan/Files",
#                 sim_config_filename="INPUT.swn",
#                 command="swanrun",
#                 on=machine_group)

# Uses swan.exe. Note: the simulation file must be called INPUT
# task = snl_swan.run(input_dir="/Path/to/My/Snl-Swan/Files",
#                 command="swan.exe",
#                 on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
