"""NWChem example."""
import inductiva

# Allocate machine
machine_group = inductiva.resources.MachineGroup("c3d-standard-180")

# Initialize the Simulator
nwchem = inductiva.simulators.NWChem()

# Run simulation with config files in the input directory
task = nwchem.run(input_dir="/path/to/my/nwchem/files",
                  sim_config_filename="my_config_file.nw",
                  on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
machine_group.terminate()

task.download_outputs()