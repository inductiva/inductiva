"""FDS example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    machine_type="c3d-standard-180",
    provider="GCP")

# Initialize the Simulator
fds = inductiva.simulators.FDS()

# Run simulation with config files in the input directory
task = fds.run(input_dir="path/to/my/fds/files",
               sim_config_filename="my_config_file.fds",
               on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
