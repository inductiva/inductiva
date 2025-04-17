"""SWAN example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highcpu-180")

# Initialize the Simulator
swan = inductiva.simulators.SWAN(version="41.45")

# Run simulation with config files in the input directory
task = swan.run( \
    input_dir="/path/to/my/swan/files",
    sim_config_filename="my_config_file.swn",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
