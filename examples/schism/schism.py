"""SCHISM example."""
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")

# Initialize the Simulator
schism = inductiva.simulators.SCHISM()

# Run simulation with config files in the input directory
task = schism.run(input_dir="/path/to/my/schism/files",
                  n_vcpus=90,
                  num_scribes=2,
                  on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
