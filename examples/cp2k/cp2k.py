""" CP2K example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
cp2k = inductiva.simulators.CP2K()

# Run simulation with config files in the input directory
task = cp2k.run(input_dir="/Path/to/My/cp2k/Files",
                sim_config_filename="my_config_file.inp",
                on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()