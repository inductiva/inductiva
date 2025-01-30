"""SPlisHSPlasH example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
splishsplash = inductiva.simulators.SplishSplash()

# Run simulation
task = splishsplash.run(input_dir="/path/to/my/splishsplash/files",
                        sim_config_filename="config.json",
                        on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
