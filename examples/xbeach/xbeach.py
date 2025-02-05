"""XBeach example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach(version="1.24")

# Run simulation with configuration files in the input directory
task = xbeach.run(input_dir="/path/to/my/xbeach/files",
                  sim_config_filename="my_config_file.txt",
                  on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()
