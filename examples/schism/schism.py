"""SCHISM example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
schism = inductiva.simulators.SCHISM()

# Run simulation with config files in the input directory
task = schism.run(input_dir="/path/to/my/schism/files",
                  num_scribes=2,
                  on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()