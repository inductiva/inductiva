"""CaNS example"""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
cans = inductiva.simulators.CaNS()

# Run simulation with config files in the input directory
task = cans.run(input_dir="path/to/my/cans/files",
                sim_config_filename="my_config_file.nml",
                on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()
